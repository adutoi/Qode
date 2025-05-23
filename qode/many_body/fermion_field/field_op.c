/*   (C) Copyright 2018, 2019, 2023 Yuhong Liu and Anthony Dutoi
 *
 *   This file is part of Qode.
 *
 *   Qode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Qode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Qode.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * The basic operational principle is that an array of configurations is given
 * along with a state vector (having the same length).  These configurations
 * are bit strings (perhaps consisting of multiple consecutive integers), where
 * each bit represents the occupation of an orbital.  The action of field
 * operators are implemented by bit flips on copies of the configuration, while
 * simultaneously keeping track of the phase (using the prepend convention and
 * the orbital ordering described below).  After action of a string of field
 * operators on a ket configuration, the location of the new configuration in
 * a list of bra configurations is looked up, in order to find the bra onto
 * which it projects.  This look-up relies on the fact that the configuration
 * strings are stored in ascending order according to the interpretation of
 * their bit strings as integers.  This step (implemented via bisection search)
 * is the most expensive part of this algorithm.
 *
 * On the inside of this code (which may differ from the outside) we imagine
 * the binary representation of the integers used to represent configurations
 * written in the usual way, with the low bit (representing the lowest-index
 * orbital) to the right.  Therefore, to make things conceptually consistent,
 * it helps to imagine arrays written right-to-left in storage (opposite the
 * natural direction for a left-to-right reader), so that the axes of the
 * integrals/density tensors associated with the same orbital sets as the bit
 * strings run in the same direction.  This is particularly convenient for the
 * *arrays* of configuration integers that need to be used to allow for systems
 * with more orbitals than can fit in a single BigInt; the lowest-order bits
 * (that is, the lowest-index orbitals) are represented in the 0-th element of
 * such an array, and so on, so it helps to think of the 1st element as being
 * to the left of the 0th, etc.  These configuration integers come in as one
 * contiguous block (to ensure they are contiguous for efficiency) and get
 * chopped up as they are used.
 *
 * A "wisdom" mechanism has been implemented that essentially allows lookup
 * tables produced in via one execution of a function to be used in a later
 * call.  This is ostensibly to speed it up, though the gains on a single
 * thread are modest and bus traffic causes it to actually slow down on mutiple
 * threads.  We nevertheless keep it around because the object itself has other
 * uses as a compact representation of a requested transition density tensor
 * between *all* pairs of configurations in the basis.  More detail is given
 * where relevant, but it is worthwhile to mention here that, in order to store
 * both index and phase of a non-zero element without wasting memory or
 * making the code even more horrible, negative indices are used to indicate
 * that an element at a certain index is negative in value.  This necessitates
 * using fortran-style indexing for the storage, starting at 1.  Zero is then
 * used as a signal that a certain operator string on a ket overlaps no bra ket
 * (it otherwise overlaps maximally one bra, whose absolute index is stored).
 *
 * Note: the typing (BigInt, PyInt, etc) has been thought through more
 * carefully than it may appear so proceed with caution.  Some functions must
 * be typed for communication with python.  Everything is ok as long as BigInt
 * is not smaller than PyInt, and both are not smaller than int.
 *
 * Ideas for making this better:
 * - allow input about what ranges of creation/annihilation have nonzero elements
 * - feed in integrals that account for frozen core while working only with active-space configs
 * - spin symmetry?
 */

#include <stdlib.h>       // exit()
#include <string.h>       // memcpy()
#include <math.h>         // fabs()
#include "PyC_types.h"    // PyInt, BigInt, Double
#ifdef _OPENMP
#include <omp.h>          // multithreading
#endif



// The operation of the recursive kernel for performing the action of an operator on a vector
// to produce a new vector or to build transition-density tensors between vectors is largely the
// same, so there is one piece of code, where one of three modes is activated upon bottoming out.
//
#define OP_ACTION   1    // for acting an operator on a vector to produce a new vector
#define COMPUTE_D   2    // for computing density tensors between sets of bras and kets
#define WISDOM_ONLY 3    // for computing the lookup tables only (only GENERATE below makes sense with this)

// What to do with the pointers to the wisdom (lookup) tables provided at the top level.
#define IGNORE      0    // ignore them (probably NULL)
#define GENERATE    1    // populate them
#define APPLY       2    // use the information to avoid calling the bisection_search function

// This utility divides val into its sign and absolute value without having to worry about
// whether to hardcode the use of abs(), labs(), llabs(), etc.  It also shifts the absolute
// value from fortran-style indexing to C-style indexing (see comments above and below).
//
void extract(int* sign, BigInt* absval, BigInt val)
    {
    if (val >= 0)
        {
        *sign   = 1;
        *absval = val - 1;
        }
    else
        {
        *sign   = -1;
        *absval = (-val) - 1;
        }
    return;
    }

// The number of orbitals represented by a single BigInt.
// Self-standing function to be accessible to python for packing the configs arrays.
//
PyInt orbs_per_configint()
    {
    return 8*sizeof(BigInt) - 1;    // pretty safe to assume a byte is 8 bits.  Less 1 bit for sign.
    }



// The index of a config in an array of configs (given that each configuration requires
// n_configint BigInts).  Returns -1 if config not in configs.  The last two arguments
// are the inclusive initial bounds.
// Also needs to be accessible to python.
//
PyInt bisect_search(BigInt* config, BigInt* configs, PyInt n_configint, PyInt lower, PyInt upper)
    {
    PyInt   i, half;        // Get used repeatedly, ...
    BigInt  deviation;      // ... so just declare once ...
    BigInt* test_config;    // ... in the old-fashioned way.

    // Note that we will only ever care if deviation is >, ==, or < zero.  Therefore,
    // we test the high-order components of two configurations first and we only test
    // the ones below it if those are equal.  Otherwise, we already have our answer.

    // Make sure config is not below the lower bound
    deviation = 0;
    test_config = configs + (lower * n_configint);    // pointer arithmetic for start of test_config
    for (i=n_configint-1; i>=0 && deviation==0; i--)  {deviation = config[i] - test_config[i];}
    if (deviation < 0)  {return -1;}

    // Make sure config is not above the upper bound
    deviation = 0;
    test_config = configs + (upper * n_configint);    // pointer arithmetic for start of test_config
    for (i=n_configint-1; i>=0 && deviation==0; i--)  {deviation = config[i] - test_config[i];}
    if (deviation > 0)  {return -1;}

    // Narrow it down to a choice of one or two numbers.  So any time that the contents of this
    // loop runs, there are at least three choices, so that half will not be one of the boundaries.
    // It can happen that this loop ends with upper==lower, which just makes the two code blocks
    // that follow redundant (but that will not lead to a wrong result).
    while (upper-lower > 1)
        {
        half = (lower + upper) / 2;
        deviation = 0;
        test_config = configs + (half * n_configint);    // pointer arithmetic for start of test_config
        for (i=n_configint-1; i>=0 && deviation==0; i--)  {deviation = config[i] - test_config[i];}
        if      (deviation == 0)  {return half;}         // only happens if all components equal
        else if (deviation  > 0)  {lower = half+1;}
        else                      {upper = half-1;}
        }

    // If the config is the lower of the two (or one) choices, return that index.
    deviation = 0;
    test_config = configs + (lower * n_configint);    // pointer arithmetic for start of test_config
    for (i=n_configint-1; i>=0 && deviation==0; i--)  {deviation = config[i] - test_config[i];}
    if (deviation == 0)  {return lower;}              // only happens if all components equal

    // If the config is the upper of the two (or one) choices, return that index.
    deviation = 0;
    test_config = configs + (upper * n_configint);    // pointer arithmetic for start of test_config
    for (i=n_configint-1; i>=0 && deviation==0; i--)  {deviation = config[i] - test_config[i];}
    if (deviation == 0)  {return upper;}              // only happens if all components equal

    // Else, the default is to announce that the config is not in the configs array.
    return -1;
    }



// The recursive kernel for looping over orbital indices of a string of field operators for the
// purpose of either acting an operator that is a linear combination of such strings or
// computing the separate matrix elements of all such strings.  See comments with the driver
// functions below for more context.
//
// The general principle of operation is the that driver code loops over the ket configurations.
// These configurations may be components of multiple ket states, but since resolving action of
// field-operator strings upon them is the most expensive thing, the loop over these configurations
// is pulled to the outside (in the driver code); this is followed by the loops over the orbital
// indices (done by recursion here), after which the action is completely resolved (done in layers
// to keep most effort out of the more inner loops); finally the loops over the ket (and perhaps
// bra) coefficents are most deeply nested (performed here when the recursion bottoms out).  The
// latter, lower orbital indices (for operators closest to the ket) are looped over at the
// outer-most level of recursion.  To increae efficiency, this code dynamically keeps track of
// which orbitals are occupied and empty as the configuration is modified with each level of
// recursion.  The tricky part, which leads to more lines of code than one might first think to
// write, is to get the recursion to bottom out at 1 loop instead of 0.  This is because things
// can be made more efficient if it is known when one is at the inner-most loop, and it is the
// inner-most loop that actually matters.
//
// The wisdom (lookup) tables, when requested are filled up on a per-ket basis.  With the orbital-
// index combinations linearly enumerated exactly as they are looped over here, the index of the
// bra onto which the action of the string onto the ket would project is stored.  The index is
// stored fortran-style (shifted up by one, as explained above), and made negative if the overlap
// is negative.  Since it is of interest to outside users, more is said about the format of this
// table (which must be accompanied by ket occupancy info) is written below with the function
// determinant_densities(), whose job it is to compute and return only these tables.
//
void resolve_recur(int      mode,              // OP_ACTION or COMPUTE_D (determines whether using Psi_L for storing new states or as bras)
                   PyInt    n_create,          // number of creation operators at present level of recursion
                   PyInt    n_annihil,         // number of annihilation operators at present level of recursion
                   Double** Psi_L,             // states being produced (LHS of equation) for OP_ACTION; states in the bra (on left) for COMPUTE_D
                   PyInt    n_Psi_L,           // number of states in Psi_L (for OP_ACTION, must have n_Psi_L==n_Psi_R, below)
                   BigInt*  configs_L,         // configuration strings representing the basis for the states in Psi_L
                   PyInt    n_configs_L,       // number of configurations in the basis configs_L
                   PyInt    n_configint_L,     // number of BigInts needed to store a single configuration in configs_L
                   Double** Psi_R,             // states being acted on (RHS of equation) for OP_ACTION; states in the ket (on right) for COMPUTE_D
                   PyInt    n_Psi_R,           // number of states in Psi_R (for OP_ACTION, must have n_Psi_L==n_Psi_R, above)
                   BigInt*  config_R,          // ket (right-hand) configuration being acted upon at present layer of recursion
                   PyInt    config_idx_R,      // index of the configuration (in the right-hand basis) acted upon at the *top* layer of recursion
                   PyInt    n_configint_R,     // number of BigInts needed to store the ket configuration being acted on, at any level of recursion
                   Double** tensors,           // tensor of matrix elements (sole entry) for OP_ACTION, or storage for output (array of arrays) for COMPUTE_D
                   PyInt    n_orbs,            // edge dimension of the tensor(s)
                   PyInt    global_phase,      // a global phase to be applied to the operator action
                   int*     occupied,          // indices of orbitals that are occupied in the configuration at the present level of recursion (not necessarily in order)
                   int      n_occ,             // number of orbitals that are occupied at the present level of recursion
                   int*     empty,             // indices of orbitals that are empty in the configuration at the present level of recursion (not necessarily in order)
                   int      n_emt,             // number of orbitals that are empty at the present level of recursion
                   int*     cum_occ,           // cumulative number of orbitals at or below a given index that are occupied at present level of recursion
                   int      permute,           // number of permutations performed so far to satisfy the field-operator prepend convention for ket orbitals in descending order
                   int      op_idx,            // recursive build of index for the tensors array
                   int      stride,            // stride to be applied at this level of recursion in order to build op_idx (must start as 1)
                   int      factor,            // recursive build of factor to avoid looping over redundant matrix elements (must start as 1)
                   int      p_0,               // initial orbital index at this level, to avoid looping over redundant matrix elements (must start as 0)
                   Double   thresh,            // perform no further work if result will be smaller than this
                   PyInt    wisdom,            // IGNORE or GENERATE (determines whether we will generate the wisdom/lookup tables, APPLY handled by other code)
                   BigInt*  wisdom_det_idx,    // for the given ket config, a list of the (possibly negated) index that each respective field-operator string gives projection onto
                   BigInt*  wisdom_op_idx)     // running index for the entries in wisdom_det_idx (done by reference because should not be reset in outer loops)
    {
    // Some admin that needs to be done for either mode at any level of recursion
    int    n_bits = orbs_per_configint();                        // number of bits/orbitals in a BigInt
    int    n_bytes_config_R = n_configint_R * sizeof(BigInt);    // number of bytes in the ket config (for memcpy)
    BigInt p_config_R[n_configint_R];                            // a place to store modified configurations
    int    p_n;                                                  // see ...
    int*   orb_list;                                             // ... below
    if (n_annihil > 0)    // will loop over annihilation (occupied) index
        {
        p_n      = n_occ;       // upper limit of the orbital loop if doing annihlation operator
        orb_list = occupied;    // resolution of counting index to orbital index draws from occupied orbitals of present configurations
        }
    else                  // will loop over creation (empty) index
        {
        p_n      = n_emt;       // upper limit of the orbital loop if doing creation operator
        orb_list = empty;       // resolution of counting index to orbital index draws from empty orbitals of present configurations
        }

    if (n_annihil + n_create > 1)    // recursive part (there is still >1 loop to go)
        {
        int  n_bytes_cum_occ = n_orbs * sizeof(int);    // number of bytes in cum_occ array (for memcpy)
        int  p_cum_occ[n_orbs];                         // a place to store modified cum_occ arrays
        int  reset_p_0 = 0;                             // see ...
        int  occ_change;                                //   ...
        int* other_orb_list_entry;                      // ... below
        if (n_annihil > 0)    // will loop over annihilation (occupied) index
            {
            factor *= n_annihil;                          // increase the redundancy factor (only used for OP_ACTION)
            n_annihil--;                                  // there will be one less annihilation loop
            other_orb_list_entry = empty + n_emt++;       // we will eventually add the orbital index to the end of the empty array (whose length is incremented)
            occ_change = -1;                              // occupancies in cum_occ array will go down
            if (n_annihil == 0) {reset_p_0 = 1;}          // if this is the last annihilation operator we will reset the beginning index of the orbital loops (for creation loops)
            }
	else                  // will loop over creation (empty) index
            {
            factor *= n_create;                           // increase the redundancy factor (only used for OP_ACTION)
            n_create--;                                   // there will be one less creation loop
            other_orb_list_entry = occupied + n_occ++;    // we will eventually add the orbital index to the end of the occupied array (whose length is incremented). Ignored (see below)
            occ_change = +1;                              // occupancies in cum_occ array will go up
            }
        for (int p_=p_0; p_<p_n; p_++)    // loop over the "counting" index for either occupieds or empties (starting from the given value; see below)
            {
            int p = orb_list[p_];    // "absolute" index of the orbital
            int Q = p / n_bits;      // Q=quotient:  in which component of config is orbital p?
            int r = p % n_bits;      // r=remainder: which bit in ^this component is this orbital?
            int p_permute = permute + cum_occ[n_orbs-1] - cum_occ[p];     // how many permutations does it take to get to/from position p to the front (prepend convention)
            memcpy(p_config_R, config_R, n_bytes_config_R);               // a copy of the original configuration ...
            p_config_R[Q] = p_config_R[Q] ^ ((BigInt)1<<r);               // ... with occupancy of postion p flipped
            memcpy(p_cum_occ, cum_occ, n_bytes_cum_occ);                  // a copy of the cum_occ array ...
            for (int i=p; i<n_orbs; i++) {p_cum_occ[i] -= occ_change;}    // ... with occupancies appropriately altered
            int q_0 = p_ + 1;             // the beginning of the next loop starts above the current index ...
            if (reset_p_0)  {q_0 = 0;}    // ... unless we are switching from annihilation to creation operators
            *other_orb_list_entry = p;    // if we annihlated orbital p, we will want to loop over its creation as well (vice versa has no effect (or harm))
            // recur, passing through appropriately modified quantities (see below about inline updates)
            resolve_recur(mode, n_create, n_annihil, Psi_L, n_Psi_L, configs_L, n_configs_L, n_configint_L, Psi_R, n_Psi_R, p_config_R, config_idx_R, n_configint_R, tensors, n_orbs, global_phase, occupied, n_occ, empty, n_emt, p_cum_occ, p_permute, op_idx+p*stride, stride*n_orbs, factor, q_0, thresh, wisdom, wisdom_det_idx, wisdom_op_idx);
            }
        }
    else if (mode == OP_ACTION)    // bottom out option
        {
        for (int p_=p_0; p_<p_n; p_++)    // final orbital loop (see above)
            {
            int p = orb_list[p_];                                   // absolute index (see above)
            Double val = factor * tensors[0][p*stride + op_idx];    // finish building tensor index (done inline with recursion above) and get integral from only tensor
            int compute_it = (fabs(val) > thresh);                  // do nothing if the integral is too small (thresh considers also ket coefficient) ...
            if (compute_it || (wisdom == GENERATE))                 // ... unless we are generating the lookup tables (then we need to at least compute the index)
                {
                int Q = p / n_bits;                                // build ...
                int r = p % n_bits;                                // ... modified configuration ...
                memcpy(p_config_R, config_R, n_bytes_config_R);    // ... as discussed ...
                p_config_R[Q] = p_config_R[Q] ^ ((BigInt)1<<r);    // ... above
                BigInt config_idx_L = bisect_search(p_config_R, configs_L, n_configint_L, 0, n_configs_L-1);    // EXPENSIVE! -- find left-basis index of full string on right config
                int i = (*wisdom_op_idx)++;                                         // get and increment the running index of the lookup table (whether used or not)
                if (wisdom == GENERATE)  {wisdom_det_idx[i] = config_idx_L + 1;}    // if generating lookup table, store the left/bra index (fortran-style indexing)
                if (config_idx_L != -1)    // do nothing if action takes outside of space of configurations
                    {
                    int p_permute = permute + cum_occ[n_orbs-1] - cum_occ[p];    // final permutation and ...
                    int phase = (p_permute%2) ? -1 : 1;                          // ... multiplication by resulting phase ...
                    if (compute_it)
                        {
                        val *= global_phase * phase;                             // ... delayed until we know operation was nonzero
                        for (int v=0; v<n_Psi_R; v++)    // for each vector in the input set ...
                            {
                            Double update = val * Psi_R[v][config_idx_R];    // ... connect the input configuration ...
                            #pragma omp atomic                               // ... (in a thread-safe way) ...
                            Psi_L[v][config_idx_L] += update;                // ... to the slot of the output
                            }
                        }
                    if (wisdom == GENERATE)  {wisdom_det_idx[i] *= phase;}       // if generating lookup table, store the phase as the sign of the index (explains fortran-style indexing)
                    }
                }
            }
        }
    else if (mode == COMPUTE_D)    // bottom out option
        {
        for (int p_=p_0; p_<p_n; p_++)    // final orbital loop (see above)
            {
            int p = orb_list[p_];                              // absolute index (see above)
            int Q = p / n_bits;                                // build ...
            int r = p % n_bits;                                // ... modified configuration ...
            memcpy(p_config_R, config_R, n_bytes_config_R);    // ... as discussed ...
            p_config_R[Q] = p_config_R[Q] ^ ((BigInt)1<<r);    // ... above
            BigInt config_idx_L = bisect_search(p_config_R, configs_L, n_configint_L, 0, n_configs_L-1);    // EXPENSIVE! -- find left-basis index of full string on right config
            int i = (*wisdom_op_idx)++;                                         // get and increment the running index of the lookup table (whether used or not)
            if (wisdom == GENERATE)  {wisdom_det_idx[i] = config_idx_L + 1;}    // if generating lookup table, store the left/bra index (fortran-style indexing)
            if (config_idx_L != -1)    // do nothing if action takes outside of space of configurations
                {
                int p_op_idx = p*stride + op_idx;                            // finish building tensor index (done inline with recursion above)
                int p_permute = permute + cum_occ[n_orbs-1] - cum_occ[p];    // final permutation and ...
                int phase = (p_permute%2) ? -1 : 1;                          // ... computation of resulting phase
                int braket = 0;                                              // initialize a running index for the bra-ket pairs
                for (int vL=0; vL<n_Psi_L; vL++)    // loop over the bra states ...
                    {
                    Double coeff_L = global_phase * phase * Psi_L[vL][config_idx_L];    // ... and get the phased coefficient of the left configuration for each state
                    if (fabs(coeff_L) > thresh)    // do nothing if left coefficient is too small (thresh considers also ket coefficient)
                        {
                        for (int vR=0; vR<n_Psi_R; vR++)    // loop over the ket states ...
                            {
                            Double update = coeff_L * Psi_R[vR][config_idx_R];    // ... and add the (phased) product of the left and right coefficients ...
                            #pragma omp atomic                                    // ... (in a thread-safe way) ...
                            tensors[braket++][p_op_idx] += update;                // ... to the precomputed location in the tensor for the corresponding bra-ket pair (which is incremented)
                            }
                        }
                    else
                        {
                        braket += n_Psi_R;                                        // make sure to increment the running index even if ket loop skipped
                        }
                    }
                if (wisdom == GENERATE)  {wisdom_det_idx[i] *= phase;}          // if generating lookup table, store the left/bra index (fortran-style indexing)
                }
            }
        }
    else if (mode == WISDOM_ONLY)    // bottom out option
        {
        for (int p_=p_0; p_<p_n; p_++)    // final orbital loop (see above)
            {
            int p = orb_list[p_];                              // absolute index (see above)
            int Q = p / n_bits;                                // build ...
            int r = p % n_bits;                                // ... modified configuration ...
            memcpy(p_config_R, config_R, n_bytes_config_R);    // ... as discussed ...
            p_config_R[Q] = p_config_R[Q] ^ ((BigInt)1<<r);    // ... above
            BigInt config_idx_L = bisect_search(p_config_R, configs_L, n_configint_L, 0, n_configs_L-1);    // EXPENSIVE! -- find left-basis index of full string on right config
            int i = (*wisdom_op_idx)++;              // get and increment the running index of the lookup table
            wisdom_det_idx[i] = config_idx_L + 1;    // store the left/bra index (fortran-style indexing)
            if (config_idx_L != -1)    // do nothing if action takes outside of space of configurations
                {
                int p_op_idx = p*stride + op_idx;                            // finish building tensor index (done inline with recursion above)
                int p_permute = permute + cum_occ[n_orbs-1] - cum_occ[p];    // final permutation and ...
                int phase = (p_permute%2) ? -1 : 1;                          // ... computation of resulting phase
                wisdom_det_idx[i] *= phase;          // store the left/bra index (fortran-style indexing)
                }
            }
        }
    else
       {
       exit(EXIT_FAILURE);    // unlikely, but just in case, an exit here will be easier to debug
       }

    return;
    }



// This performs the same basic function as resolve_recur(), but uses the wisdom/lookup tables rather than a bisection search to find
// the indices of the bra states projected onto.  The APPLY value of the wisdom argument in calling code directs to here.
// There is an unfortunate amount of redundant code here but it has to be this way because we want to do the same thing faster (skipping steps).
//
void resolve_recur_wise(int      mode,              // OP_ACTION or COMPUTE_D (determines whether using Psi_L for storing new states or as bras)
                        PyInt    n_create,          // number of creation operators at present level of recursion
                        PyInt    n_annihil,         // number of annihilation operators at present level of recursion
                        Double** Psi_L,             // states being produced (LHS of equation) for OP_ACTION; states in the bra (on left) for COMPUTE_D
                        PyInt    n_Psi_L,           // number of states in Psi_L (for OP_ACTION, must have n_Psi_L==n_Psi_R, below)
                        Double** Psi_R,             // states being acted on (RHS of equation) for OP_ACTION; states in the ket (on right) for COMPUTE_D
                        PyInt    n_Psi_R,           // number of states in Psi_R (for OP_ACTION, must have n_Psi_L==n_Psi_R, above)
                        PyInt    config_idx_R,      // index of the configuration (in the right-hand basis) acted upon at the *top* layer of recursion
                        Double** tensors,           // tensor of matrix elements (sole entry) for OP_ACTION, or storage for output (array of arrays) for COMPUTE_D
                        PyInt    n_orbs,            // edge dimension of the tensor(s)
                        PyInt    global_phase,      // a global phase to be applied to the operator action
                        int*     occupied,          // indices of orbitals that are occupied in the configuration at the present level of recursion (not necessarily in order)
                        int      n_occ,             // number of orbitals that are occupied at the present level of recursion
                        int*     empty,             // indices of orbitals that are empty in the configuration at the present level of recursion (not necessarily in order)
                        int      n_emt,             // number of orbitals that are empty at the present level of recursion
                        int      op_idx,            // recursive build of index for the tensors array
                        int      stride,            // stride to be applied at this level of recursion in order to build op_idx (must start as 1)
                        int      factor,            // recursive build of factor to avoid looping over redundant matrix elements (must start as 1)
                        int      p_0,               // initial orbital index at this level, to avoid looping over redundant matrix elements (must start as 0)
                        Double   thresh,            // perform no further work if result will be smaller than this
                        BigInt*  wisdom_det_idx,    // for the given ket config, a list of the (possibly negated) index that each respective field-operator string gives projection onto
                        BigInt*  wisdom_op_idx)     // running index for the entries in wisdom_det_idx (done by reference because should not be reset in outer loops)
    {
    // Some admin that needs to be done for either mode at any level of recursion
    int    p_n;                                                  // see ...
    int*   orb_list;                                             // ... below
    if (n_annihil > 0)    // will loop over annihilation (occupied) index
        {
        p_n      = n_occ;       // upper limit of the orbital loop if doing annihlation operator
        orb_list = occupied;    // resolution of counting index to orbital index draws from occupied orbitals of present configurations
        }
    else                  // will loop over creation (empty) index
        {
        p_n      = n_emt;       // upper limit of the orbital loop if doing creation operator
        orb_list = empty;       // resolution of counting index to orbital index draws from empty orbitals of present configurations
        }

    if (n_annihil + n_create > 1)    // recursive part (there is still >1 loop to go)
        {
        int  reset_p_0 = 0;                             // see ...
        int* other_orb_list_entry;                      // ... below
        if (n_annihil > 0)    // will loop over annihilation (occupied) index
            {
            factor *= n_annihil;                          // increase the redundancy factor (only used for OP_ACTION)
            n_annihil--;                                  // there will be one less annihilation loop
            other_orb_list_entry = empty + n_emt++;       // we will eventually add the orbital index to the end of the empty array (whose length is incremented)
            if (n_annihil == 0) {reset_p_0 = 1;}          // if this is the last annihilation operator we will reset the beginning index of the orbital loops (for creation loops)
            }
	else                  // will loop over creation (empty) index
            {
            factor *= n_create;                           // increase the redundancy factor (only used for OP_ACTION)
            n_create--;                                   // there will be one less creation loop
            other_orb_list_entry = occupied + n_occ++;    // we will eventually add the orbital index to the end of the occupied array (whose length is incremented). Ignored (see below)
            }
        for (int p_=p_0; p_<p_n; p_++)    // loop over the "counting" index for either occupieds or empties (starting from the given value; see below)
            {
            int p = orb_list[p_];         // "absolute" index of the orbital
            int q_0 = p_ + 1;             // the beginning of the next loop starts above the current index ...
            if (reset_p_0)  {q_0 = 0;}    // ... unless we are switching from annihilation to creation operators
            *other_orb_list_entry = p;    // if we annihlated orbital p, we will want to loop over its creation as well (vice versa has no effect (or harm))
            // recur, passing through appropriately modified quantities (see below about inline updates)
            resolve_recur_wise(mode, n_create, n_annihil, Psi_L, n_Psi_L, Psi_R, n_Psi_R, config_idx_R, tensors, n_orbs, global_phase, occupied, n_occ, empty, n_emt, op_idx+p*stride, stride*n_orbs, factor, q_0, thresh, wisdom_det_idx, wisdom_op_idx);
            }
        }
    else if (mode == OP_ACTION)    // bottom out option
        {
        for (int p_=p_0; p_<p_n; p_++)    // final orbital loop (see above)
            {
            int p = orb_list[p_];                                   // absolute index (see above)
            Double val = factor * tensors[0][p*stride + op_idx];    // finish building tensor index (done inline with recursion above) and get integral from only tensor
            if (fabs(val) > thresh)    // do nothing if the integral is too small (thresh considers also ket coefficient)
                {
                int i = (*wisdom_op_idx)++;    // get and increment the running index of the lookup table
                int phase;                     // accumulated permutational phase
                BigInt config_idx_L;           // left-basis index of full string on right config
                extract(&phase, &config_idx_L, wisdom_det_idx[i]);    // extract the phase and bra index from the lookup table
                if (config_idx_L != -1)    // do nothing if action takes outside of space of configurations
                    {
                    val *= global_phase * phase;    // rephase the tensor element
                    for (int v=0; v<n_Psi_R; v++)    // for each vector in the input set ...
                        {
                        Double update = val * Psi_R[v][config_idx_R];    // ... connect the input configuration ...
                        #pragma omp atomic                               // ... (in a thread-safe way) ...
                        Psi_L[v][config_idx_L] += update;                // ... to the slot of the output
                        }
                    }
                }
            else
                {
                (*wisdom_op_idx)++;    // make sure to increment the running index even if loop skipped
                }
            }
        }
    else if (mode == COMPUTE_D)    // bottom out option
        {
        for (int p_=p_0; p_<p_n; p_++)    // final orbital loop (see above)
            {
            int p = orb_list[p_];          // absolute index (see above)
            int i = (*wisdom_op_idx)++;    // get and increment the running index of the lookup table
            int phase;                     // accumulated permutational phase
            BigInt config_idx_L;           // left-basis index of full string on right config
            extract(&phase, &config_idx_L, wisdom_det_idx[i]);    // extract the phase and bra index from the lookup table
            if (config_idx_L != -1)    // do nothing if action takes outside of space of configurations
                {
                int p_op_idx = p*stride + op_idx;                            // finish building tensor index (done inline with recursion above)
                int braket = 0;                                              // initialize a running index for the bra-ket pairs
                for (int vL=0; vL<n_Psi_L; vL++)    // loop over the bra states ...
                    {
                    Double coeff_L = phase * Psi_L[vL][config_idx_L];    // ... and get the phased coefficient of the left configuration for each state
                    if (fabs(coeff_L) > thresh)    // do nothing if left coefficient is too small (thresh considers also ket coefficient)
                        {
                        for (int vR=0; vR<n_Psi_R; vR++)    // loop over the ket states ...
                            {
                            Double update = coeff_L * Psi_R[vR][config_idx_R];    // ... and add the (phased) product of the left and right coefficients ...
                            #pragma omp atomic                                    // ... (in a thread-safe way) ...
                            tensors[braket++][p_op_idx] += update;                // ... to the precomputed location in the tensor for the corresponding bra-ket pair (which is incremented)
                            }
                        }
                    else
                        {
                        braket += n_Psi_R;                                        // make sure to increment the running index even if ket loop skipped
                        }
                    }
                }
            }
        }
    else
       {
       exit(EXIT_FAILURE);    // unlikely, but just in case, an exit here will be easier to debug
       }

    return;
    }



// The set-up/driver for the recursive kernel, which is itself called by more specialized code to
// perform one of two actions (described shortly).  See comments with the recursive kernel above
// for more information on the implementation, and see comments with the top-level functions below
// for more information on use cases.
//
// This function (via the recursive kernel) loops over all orbital indices associated with a given
// vacuum-normal-ordered string of field operators, for the purpose of either acting an operator
// that is a linear combination of such strings onto a set of states, or computing the separate
// matrix elements of all such strings (ie, a transition density tensor).  Which action is taken
// is determined by the first mode argument.  In the former case, the left-hand states are storage
// for the action of the operator on the right-hand states (they are incremented by the result),
// and the tensors input contains matrix elements that define the operator.  In the latter case,
// the left-hand states are interpreted as bra states for the transition density tensors, and the
// tensors argument provides storage for the result (again, incremented).  For efficiency, this
// only considers strings with descending-order combinations of orbital indices for the separate
// substrings of creation and annihilation operators (written left to right, with the ket on the
// right).  This assumes that the input integrals are antisymmetric or that the output densities
// will be subsequently antisymmetrized.  This driver function specifically performs the outermost
// loop over the ket configurations, with loops over orbital indices and ket (and perhaps bra)
// coefficients delegated to the recursive kernel.
//
void resolve(int      mode,               // OP_ACTION or COMPUTE_D (determines whether using Psi_L for storing new states or as bras)
             PyInt    n_create,           // number of creation operators
             PyInt    n_annihil,          // number of annihilation operators
             Double** Psi_L,              // states being produced (LHS of equation) for OP_ACTION; states in the bra (on left) for COMPUTE_D
             PyInt    n_Psi_L,            // number of states in Psi_L (for OP_ACTION, must have n_Psi_L==n_Psi_R, below)
             BigInt*  configs_L,          // configuration strings representing the basis for the states in Psi_L
             PyInt    n_configs_L,        // number of configurations in the basis configs_L
             PyInt    n_configint_L,      // number of BigInts needed to store a single configuration in configs_L
             Double** Psi_R,              // states being acted on (RHS of equation) for OP_ACTION; states in the ket (on right) for COMPUTE_D
             PyInt    n_Psi_R,            // number of states in Psi_R (for OP_ACTION, must have n_Psi_L==n_Psi_R, above)
             BigInt*  configs_R,          // configuration strings representing the basis for the states in Psi_R
             PyInt    n_configs_R,        // number of configurations in the basis configs_R
             PyInt    n_configint_R,      // number of BigInts needed to store a single configuration in configs_R
             Double** tensors,            // tensor of matrix elements (sole entry) for OP_ACTION, or storage for output (array of arrays) for COMPUTE_D
             PyInt    n_orbs,             // edge dimension of the tensor(s)
             PyInt    phase,              // a global phase to be applied to the operator action
             PyFloat  thresh,             // perform no further work if result will be smaller than this
             PyInt    wisdom,             // IGNORE, GENERATE, or APPLY (determines what, if anything, we will do with pointers to the wisdom/lookup tables)
             Int**    wisdom_occupied,    // for each ket config, a list (in ascending order) of the orbitals occupied in that ket
             BigInt** wisdom_det_idx,     // for each ket config, a list of the (possibly negated) index that each respective field-operator string gives projection onto
             PyInt    n_threads)          // number of threads to spread the work over
    {
    omp_set_num_threads(n_threads);       // declare the number of threads to use
    int n_bits = orbs_per_configint();    // number of bits/orbitals in a BigInt

    #pragma omp parallel for               // divide the outermost loop over the threads (omp atomic at update to avoid race condition)
    for (PyInt n=0; n<n_configs_R; n++)    // loop over configurations in ket basis
        {
        Double biggest = 0;    // eventual value of the biggest coefficient for a given ket configuration
        for (int v=0; v<n_Psi_R; v++)    // loop over ket states
            {
            Double size = fabs(Psi_R[v][n]);          // get the configuration coefficient ...
            if (size > biggest)  {biggest = size;}    // ... for each state and compare to others
            }

        if ((biggest > thresh) || (wisdom == GENERATE))    // do nothing if the configuration has no significant coefficients, unless we are generating wisdom/lookup tables
            {
            BigInt* config = configs_R + (n * n_configint_R);    // config[] is now an array of integers collectively holding the present configuration

            int occupied[n_orbs];   // for indices of orbitals that are occupied in the present configuration (will be appended out of order in recursion)
            int empty[n_orbs];      // for indices of orbitals that are empty    in the present configuration (will be appended out of order in recursion
            int cum_occ[n_orbs];    // for the cumulative number of orbitals at or below a given index that are occupied (for phase calculations)
            int n_occ = 0;          // eventual number of occupied orbitals (also a running index for cataloging their indices below)
            int n_emt = 0;          // eventual number of empty    orbitals (also a running index for cataloging their indices below)
            for (int p=0; p<n_orbs; p++)    // loop over all orbitals
                {
                int Q = p / n_bits;                                          // Q=quotient:  in which component of config is orbital p?
                int r = p % n_bits;                                          // r=remainder: which bit in ^this component is this orbital?
                if (config[Q] & ((BigInt)1<<r))  {occupied[n_occ++] = p;}    // if bit r is "on" in component Q, it is occupied, ...
                else                             {   empty[n_emt++] = p;}    // ... otherwise it is empty
                cum_occ[p] = n_occ;                                          // set after incrementing n_occ (so cumulative occupancy "counting this orb")
                }

            BigInt* wisdom_det_idx_n = (BigInt*)NULL;    // if wisdom generated or used, the row of the lookup table for this ket config (NULL otherwise)
            BigInt  wisdom_op_idx = 0;                   // if wisdom generated or used, the running index for the lookup table (intialized here, passed by reference for incrementing)
            if (wisdom == GENERATE)
                {
                for (int i=0; i<n_occ; i++)  {wisdom_occupied[n][i] = occupied[i];}    // "permanent" record of the occupied orbitals for this ket config (as opposed to the working version being copied)
                }
            if ((wisdom == GENERATE) || (wisdom == APPLY))
                {
                wisdom_det_idx_n = wisdom_det_idx[n];                                  // assign the correct pointer if used
                }

            // begin the recursive kernel that loops over orbital indices for the operator string, and coefficients for the ket (and perhaps bra) state(s)
            // dividing thresh/biggest yields a an effective threshold for multiplier of a ket coefficient (like a matrix element or a bra coefficient)
            if (wisdom == APPLY)    // use the wisdom/lookup tables
                {
                resolve_recur_wise(mode, n_create, n_annihil, Psi_L, n_Psi_L, Psi_R, n_Psi_R, n, tensors, n_orbs, phase, occupied, n_occ, empty, n_emt, 0, 1, 1, 0, thresh/biggest, wisdom_det_idx_n, &wisdom_op_idx);
                }
            else                    // find bra indices by modifying configurations and then searching, perhaps generating wisdom/lookup tables
                {
                resolve_recur(mode, n_create, n_annihil, Psi_L, n_Psi_L, configs_L, n_configs_L, n_configint_L, Psi_R, n_Psi_R, config, n, n_configint_R, tensors, n_orbs, phase, occupied, n_occ, empty, n_emt, cum_occ, 0, 0, 1, 1, 0, thresh/biggest, wisdom, wisdom_det_idx_n, &wisdom_op_idx);
                }
            }
        }

    return;
    }



// This function takes an n-electron (n_elec) operator, whose matrix elements are given (op) in a basis
// with a given number of orbitals (n_orb), and acts it on a set of a certain number (n_Psi) of state
// vectors (Psi) to produce (technically increment) another set of vectors (opPsi).  The vectors are
// linear combinations of a certain number (n_configs) of configurations (configs), where each
// configuration is represented by a bit-string composed of a number (n_configint) of BigInt integers
// whose bits refer to the same orbital basis (as discussed in the global comments). For efficiency, a
// threshold (thresh) can be given for neglecting a result, and the work can be spread over an number
// (n_threads) of threads.
//
// The number of axes of the tensor of matrix elements (integrals) depends on the electron-order of the
// operator (2 axes for 1 electron, 4 for 2 electrons, etc).  The format of the integrals tensor should
// correspond to an operator definition of the form (for two electrons)
//     operator  =  sum_pqrs V_pqrs c_p c_q a^r a^s
// where c_p and a^p are creation and annihilation operators, respectively, on orbital p, and V is the
// integrals tensor (with any prefactor is rolled into the definition of V).  The integrals tensor is
// presumed antisymmetrized, but technically, this means that the only elements ever accessed are those
// where the indices are in descending order from left to right separately for the creation and
// annihilation substrings.  Due to the fundamental antisymmetry of the operators, any operator with a
// non-antisymmetrized tensor can be mapped onto this and it will be more efficient.  A global phase
// (phase) can be given to adjust to the external convention of the relative ordering of tensor and
// field-operator string indices.
//
// Each state vector is contiguous in memory, but they need not be adjacent to each other.  The length
// of these vectors should be the same as the number of configurations.  The multi-integers used to store
// the configurations are packed one-after-another into a contiguous block of storage (thought of
// as being arranged right-to-left as discussed in the global comments).
//
// Lookup tables (wisdom...) to be used to simplify/speed future calls can also optionally be created
// or used.  See description with determinant_densities() function for more detail.
//
void op_Psi(PyInt    n_elec,             // electron order of the operator
            Double*  op,                 // tensor of matrix elements (integrals), assumed antisymmetrized
            PyInt    n_orbs,             // edge dimension of the integrals tensor
            PyInt    phase,              // a global phase to be applied to the operator action
            Double** opPsi,              // array of row vectors: incremented by output
            Double** Psi,                // array of row vectors: input vectors to act on
            PyInt    n_Psi,              // how many vectors we are acting on and producing simultaneously in Psi and opPsi
            BigInt*  configs,            // configuration strings representing the basis for the states in Psi and opPsi (see global comments above about format)
            PyInt    n_configs,          // number of configurations in the configs basis (call signature ok if PyInt not longer than BigInt)
            PyInt    n_configint,        // number of BigInts needed to store a single configuration in configs
            PyFloat  thresh,             // perform no further work if result will be smaller than this
            PyInt    wisdom,             // IGNORE, GENERATE, or APPLY (determines what, if anything, we will do with pointers to the wisdom/lookup tables)
            Int**    wisdom_occupied,    // for each ket config, a list (in ascending order) of the orbitals occupied in that ket
            BigInt** wisdom_det_idx,     // for each ket config, a list of the (possibly negated) index that each respective field-operator string gives projection onto
            PyInt    n_threads)          // number of threads to spread the work over
    {
    // call the generic driver in operator-action mode
    resolve(OP_ACTION, n_elec, n_elec, opPsi, n_Psi, configs, n_configs, n_configint, Psi, n_Psi, configs, n_configs, n_configint, &op, n_orbs, phase, thresh, wisdom, wisdom_occupied, wisdom_det_idx, n_threads);
    return;
    }

// This function takes the description of a field operator string in terms of the numbers (n_create
// and n_annihil) of creation and annihilation operators, normal ordered with respect to the vacuum,
// a set (bras) of a certain number (n_bras) of bra state vectors, and a set (kets) of a certain
// number (n_kets) of ket state vectors.  These vectors are linear combinations of numbers (n_configs_bra
// and n_configs_ket) of configurations (configs_bra and configs_ket), where each configuration is
// represented by a bit-string composed of a number (n_configint_bra and n_configint_ket) of BigInt
// integers whose bits refer to a basis of a number (n_orbs) of orbitals.
// In the array of allocated storage (rho), which should contain n_bras * n_kets pointers to
// non-adjacent (but internally contiguous) blocks of size n_orbs ** (n_create + n_annihil),
// elements of the transiton-density tensors of the following form are stored
//     D_pqr  =  <Psi^I| c_p c_q a^r |Psi_J>
// (for example, for bra I and ket J) for values of the orbitial indices p, q, and r that correspond
// to the orbitals of the configuration bit-string representations.  c_p and a^p are creation and
// annihilation operators, respectively, on orbital p.  These tensors are sorted into the linear array
// rho first by bra index (outer index) and then ket index (inner index).  A threshold (thresh) can be
// given for neglecting a result, and the work can be spread over an number (n_threads) of threads.
// A global phase (phase) can be given to adjust to the external convention of the ordering of the
// field-operator string indices.
//
// The transition density tensors only have values filled in where the indices are in descending order
// from left to right separately for the creation and annihilation substrings.  Due to the fundamental
// antisymmetry of the operators, the missing elements are redundant with these (to within a phase) and
// should be populated by a antisymmetrization step. This is more efficient that computing them independently.
//
// Each bra and ket state vector is contiguous in memory, but they need not be adjacent to each other.
// The length of these vectors should be the same as the number of their respective configurations.  The
// multi-integers used to store the configurations are packed one-after-another into a contiguous block of
// storage (thought of as being arranged right-to-left as discussed in the global comments).
//
// Lookup tables (wisdom...) to be used to simplify/speed future calls can also optionally be created
// or used.  See description with determinant_densities() function for more detail.
//
void densities(PyInt    n_create,           // number of creation operators
               PyInt    n_annihil,          // number of annihilation operators
               Double** rho,                // array of storage for density tensors (for each bra-ket pair)
               PyInt    n_orbs,             // edge dimension of each density tensor
               PyInt    phase,              // a global phase to be applied to the operator action
               Double** bras,               // array of row vectors: bras for transition-density tensors
               PyInt    n_bras,             // number of bras
               BigInt*  configs_bra,        // configuration strings representing the basis for the bras (see global comments above about format)
               PyInt    n_configs_bra,      // number of configurations in the basis configs_bra (call signature ok if PyInt not longer than BigInt)
               PyInt    n_configint_bra,    // number of BigInts needed to store a single configuration in configs_bra
               Double** kets,               // array of row vectors: kets for transition-density tensors
               PyInt    n_kets,             // number of kets
               BigInt*  configs_ket,        // configuration strings representing the basis for the kets (see global comments above about format)
               PyInt    n_configs_ket,      // number of configurations in the basis configs_ket (call signature ok if PyInt not longer than BigInt)
               PyInt    n_configint_ket,    // number of BigInts needed to store a single configuration in configs_ket
               PyFloat  thresh,             // perform no further work if result will be smaller than this
               PyInt    wisdom,             // IGNORE, GENERATE, or APPLY (determines what, if anything, we will do with pointers to the wisdom/lookup tables)
               Int**    wisdom_occupied,    // for each ket config, a list (in ascending order) of the orbitals occupied in that ket
               BigInt** wisdom_det_idx,     // for each ket config, a list of the (possibly negated) index that each respective field-operator string gives projection onto
               PyInt    n_threads)          // number of threads to spread the work over
    {
    // call the generic driver in compute-densities mode
    resolve(COMPUTE_D, n_create, n_annihil, bras, n_bras, configs_bra, n_configs_bra, n_configint_bra, kets, n_kets, configs_ket, n_configs_ket, n_configint_ket, rho, n_orbs, phase, thresh, wisdom, wisdom_occupied, wisdom_det_idx, n_threads);
    return;
    }

// This function produces only the wisdom/lookup tables for raw use by the user because they can
// also be interpreted as a list of all non-zero elements of a given transition density tensor
// between *all* bra and ket configurations in the basis.
//
// It takes the description of a field operator string in terms of the numbers (n_create and
// n_annihil) of creation and annihilation operators, normal ordered with respect to the vacuum,
// and sets of bra and ket configurations (configs_bra and configs_ket, of number n_configs_bra
// and n_configs_ket).  Each configuration is represented by a bit-string composed of a number
// (n_configint_bra and n_configint_ket) of BigInt integers whose bits refer to a basis of a
// number (n_orbs) of orbitals.
//
// Each ket gets one "row" in each of the lookup tables (wisdom_occupied and wisdom_det_idx).
// One row of wisdom_occupied will give (in ascending order) the indices of the orbitals that
// are occupied in a given ket configuration.  For this configuration, the corresponding row
// of wisdom_det_idx gives the index of the bra into which each combination of orbital indices
// gives nonzero projection (of either +1 or -1).  These indices are in fortran-style convention
// (starting at 1 instead of 0) and have negative value if the projection negative.  If the
// value is zero, it means that the operator string acting on the ket does not belong to the
// set of bra configurations.  The main consideration is then the ordering of the combinations
// of orbital indices.  These are best described in a looping structure
//     - outer "loops" corespond to operators closer to the ket (so, to the right, not left,
//       annihilation before creation)
//     - only combinations of where the creation and annihilation substrings are each in
//       decending order (left to right, so ascending in the loop structure) are considered
//       since antisymmetry can be handled externally
//     - for the creation operators, after an aforementioned ascending-order loop over orbitals
//       originally empty in the ket is complete, the indices of previously annihilated orbitals
//       are appended, in the order they were annihilated (lower indices closer to the ket were
//       annihilated first and will created first, again closer to the ket).
//
// The work can be spread over an number (n_threads) of threads.  The multi-integers used to store
// the configurations are packed one-after-another into a contiguous block of storage (thought of
// as being arranged right-to-left as discussed in the global comments).  The use of the (global)
// phase argument above is to make sure that the population of this is valid as wisdom in other
// functions, while still allowing for convention adjustments with the outside world.
//
void generate_wisdom(PyInt    n_create,           // number of creation operators
                     PyInt    n_annihil,          // number of annihilation operators
                     PyInt    n_orbs,             // edge dimension of each density tensor
                     BigInt*  configs_bra,        // configuration strings representing the basis for the bras (see global comments above about format)
                     PyInt    n_configs_bra,      // number of configurations in the basis configs_bra (call signature ok if PyInt not longer than BigInt)
                     PyInt    n_configint_bra,    // number of BigInts needed to store a single configuration in configs_bra
                     BigInt*  configs_ket,        // configuration strings representing the basis for the kets (see global comments above about format)
                     PyInt    n_configs_ket,      // number of configurations in the basis configs_ket (call signature ok if PyInt not longer than BigInt)
                     PyInt    n_configint_ket,    // number of BigInts needed to store a single configuration in configs_ket
                     Int**    wisdom_occupied,    // for each ket config, a list (in ascending order) of the orbitals occupied in that ket
                     BigInt** wisdom_det_idx,     // for each ket config, a list of the (possibly negated) index that each respective field-operator string gives projection onto
                     PyInt    n_threads)          // number of threads to spread the work over
    {
    // call the generic driver in compute-densities mode
    resolve(WISDOM_ONLY, n_create, n_annihil, (Double**)NULL, 0, configs_bra, n_configs_bra, n_configint_bra, (Double**)NULL, 0, configs_ket, n_configs_ket, n_configint_ket, (Double**)NULL, n_orbs, 1, 0., GENERATE, wisdom_occupied, wisdom_det_idx, n_threads);
    return;
    }





/*
// could come in handy again some day . . .
#include <stdio.h>        // printf() for debug
void print_config(BigInt* config, PyInt n_orbs, int n_bits, char* beg, char* end)
    {
    printf(beg);
    for (int p=n_orbs-1; p>=0; p--)    // loop over all orbitals
        {
        int Q = p / n_bits;                                // Q=quotient:  in which component of config is orbital p?
        int r = p % n_bits;                                // r=remainder: which bit in ^this component is this orbital?
        if (config[Q] & ((BigInt)1<<r))  {printf("1");}    // if bit r is "on" in component Q, it is occupied, ...
        else                             {printf("0");}    // ... otherwise it is empty
        }
    printf(end);
    return;
    }
*/
