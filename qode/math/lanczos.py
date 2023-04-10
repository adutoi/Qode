#    (C) Copyright 2018 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#
import time
import numpy
from .space      import conjugate, sqrt, linear_inner_product_space
from .vector_set import vector_set

# Warning, this has only been dubugged for the real-symmetric case.



def _process_new_vec(v,vX,b,reorthonormalize):
    """\
    This function finishes one iteration in building the following matrix
    projection of a Hermitian linear operator H, per the band-Lanczos algorithm,
    which presumes that the procedure is initialized with a set of B orthonormal
    vectors.

             .             |-2B>      |-(B+n)>   |-(B+1)>   |-B>        |-n>        |-1>       |X>
                  .
                       .
    <-2B|                  b[B][-B]   b[n][-B]*  b[1][-B]*  b[0][-B]    0           0          0
    <-(B+n)|               b[n][-B]   b[B][n]    b[n][ n]*  b[1][ n]*   b[0][ n]    0          0
    <-(B+1)|               b[1][-B]   b[n][n]    b[B][-1]   b[n][-1]*   b[1][-1]*   b[0][-1]   0
    <-B|                   b[0][-B]   b[1][n]    b[n][-1]   b_BX        b*_nX       b*_1X      b_0X
    <-n|                   0          b[0][n]    b[1][-1]   b_nX        -           -          -
    <-1|                   0          0          b[0][-1]   b_1X        -           -          -
    <X|                    0          0          0          b_0X        -           -          -

    The argument b is a list of B+1 lists, where list b[B] is the diagonal of the
    matrix being built and b[0] is the lowest subdiagonal.  On entry, the [-1] element
    of each sublist is one of the most recently computed column of subdiagonals.

    On entry, the v[-1] (or |-1>) is the most recently computed addition to the Lanczos set
    that spans the direct sum of the Krylov spaces.  In this notation, on entry, vX is |X~~> = H|-B>,
    which is the vector that undergoes processing here to become |X>, the next addition
    to the Lanczos set.  H is a linear operator, whose action is completed before this is
    called.  The finished orthonormalized vector is then appended to v in place, as are the
    newly computed matrix elements b_nX.

    The bottom right matrix block is still undefined at at the end of the
    iteration.  It is left up to the calling routine whether any of the old vectors
    are released.

    As this function call is a "primitive", there are som non-intuitive expectations on
    the arguments.  Vectors |-2B> through |-1> must be present for all iterations.
    On the first iteration, it is expected that the latter B vectors are the starting
    block and that the former B are all zero.  Similarly, b[0] should be initialized
    as a list of B zeros, b[1] should be initialized with B-1 zeros, and so forth,
    where b[B] begins as an empty list.

    If the reorthonormalize option is used, then the Lanczos vector created in this
    iteration is Gram-Schmidt orthogonalized to all existing Lanczos vectors
    (and renormalized); this is redundant, in principle, but it is useful if |-B>
    is nearly an eigenvector as input, which can occur near convergence of some algorithms.
    """
    debug = False
    B = len(b) - 1	# block size
    # compute/gather matrix elements b_nX = <-n|H|-B> = <-n|X~~>
    b_nX = [0]						# Placeholder for b_0X
    for n in range(1,B+1):  b_nX += [(v[-n]|vX).real]	# Should I check before casting?  H is supposed to be Hermitian ...
    for n in range(1,B+1):  b_nX += [conjugate(b[B-n][-n])]
    # orthogonalize:  |X~>   =   |X~~> - sum_{n=1 to 2B} |-n><-n|X~~>   =   |X~~> - sum_{n=1 to 2B} |-n>b_nX
    for n in range(1,2*B+1):
        if debug:  print("|X~~> -= {} |-{}>".format(b_nX[n],n))
        vX -= b_nX[n] * v[-n]
    # compute final subdiagonal:  b_0X = <X|H|-B> = sqrt(<X~|X~>)
    b_nX[0] = sqrt(vX|vX)
    if debug:  print("Old subdiagonals and new additions:\n",b,"\n",b_nX)
    # normalize:  |X> = |X~> / sqrt(<X~|X~>)
    vX /= b_nX[0]
    # optional Gram-Schmidt
    if reorthonormalize:
        for vI in v:  vX -= vI * (vI|vX)
        vX /= sqrt(vX|vX)
    # update lists
    for n in range(B+1):  b[n] += [b_nX[n]]
    v += [vX]

def _iteration(H,v,b,reorthonormalize):
    """\
    This function performs one iteration in building the matrix described in
    _process_new_vec.  It performs the action of a linear operator H onto the
    vector v[-B], where B is the block size obtained from the dimension of b.
    This can be denoted |X~~> = H|-B>.  The vector |X~~> is processed into 
    the next Lanczos vector |X>, meanwhile building the corresponding components
    of the Lanczos projection of H, as descrived above.
    """
    debug = False
    B = len(b) - 1	# block size
    if debug:
        print("Entering lanczos._iteration with B={}, reorthonormalize={}.".format(B,reorthonormalize))
        print("Vector addresses:\n",list(v))
    # create:  |X~~> = H|-B>
    vX = H | v[-B]					# This is the most expensive step.
    _process_new_vec(v,vX,b,reorthonormalize)
    if debug:
        print("New subdiagonals:\n",b)
        print("New vector addresses:\n",list(v))
        print("Exiting lanczos._iteration")

def _block_iteration(H,v,b,n,reorthonormalize):
    """\
    This function performs n iterations in building the matrix described in
    _process_new_vec.  It performs the action of a linear operator H onto the
    vector v[-B], and the n-1 vectors that follow it, up to v[-B+n-1], where B
    is the block size obtained from the dimension of b.  These can be denoted
    as successive values of |X~~>, which are processed into the next Lanczos vectors |X>,
    meanwhile building the corresponding components of the Lanczos projection of H, as descrived above.
    """
    debug = False
    B = len(b) - 1	# block size
    if n>B:  raise Exception("block action cannot exceed size of initial Lanczos block")	# should never happen, since only called from other functions in here
    if debug:
        print("Entering lanczos._block_iteration with B={}, n={}, reorthonormalize={}.".format(B,n,reorthonormalize))
        print("Vector addresses:\n",list(v))
    # create:  |X~~> = H|Y> for Y in (-B ... -B+n-1)
    vX = H.act_on_vec_block(v[-B:(-B+n or None)])						# Someday wan this to be    vX = H | v[-B:(-B+n or None)] ... for that, space needs to understand vector_sets as well as vectors
    for vXi in vX:  _process_new_vec(v,vXi,b,reorthonormalize)
    if debug:
        print("New subdiagonals:\n",b)
        print("New vector addresses:\n",list(v))
        print("Exiting lanczos._block_iteration")



def projection(H,v,n,block_action=None,save=True,autocomplete=False,reorthonormalize=True):
    """\
    This function builds and returns a square matrix projection of the linear operator H
    as a numpy matrix of dimension n x n, given a starting block of B orthonormal vectors v.

    Orthonormality of the input v is enforced by an in place Gram-Schmidt of the vectors therein.
    If block_action is set to an integer, it creates the Lanczos vectors by acting on block_action
    number of vectors at a time (or less in the case of a nonzero remainder to n/block_action ... 
    block_action cannot be greater than B, and if it is, it is quietly set equal to B at the beginning
    of the algorithm).

    If the save option is True, then a tuple is returned in which the first member is the projection,
    and the second is a list of vectors generated.  If the autocomplete option is used, then a projection
    of dimension B larger than known explicitly is returned, where the bottom-right matrix block
    is approximated using the diagonals from the block that precedes it.  If the reorthonormalize is True,
    then each Lanczos vector is Gram-Schmidt orthogonalized to all currently stored Lanczos vectors
    (and renormalized) at creation time; this is redundant, in principle, but it is useful if any of
    the vectors in v is nearly an eigenvector, which can occur near convergence of some algorithms.
    that use this.
    """
    debug = False		# Set to true to turn on verbose printing
    v = list(v)			# might come in as any iterable ... also do not want to touch original list (though contents may change due to initial Gram-Schmidt)
    B = len(v)			# block size
    if (block_action is not None) and (block_action>B):  block_action = B	# make sure we never try to act on more vectors than allowed
    if debug:  print("Entering lanczos.projection with B={}, n={}, block_action={}, save={}, autocomplete={}, and reorthonormalize={}.".format(B,n,block_action,save,autocomplete,reorthonormalize))
    # Enforce orthonormality of input block
    for i,Vi in enumerate(v):
        for Vj in v[:i]:  Vi -= Vj * (Vj|Vi)
        Vi /= sqrt(Vi|Vi)
    # initialize arrays for diagonals and subdiagonals (b) and Lanczos vectors (v)
    b = [[]]
    zeros = []
    for _ in range(B):
        zeros += [0]
        b = [list(zeros)] + b		# "Extra" elements here are expected to exist in a "dumb" iteration (and need to copy list) ...
    v = vector_set(H.space, zeros+v)	#  ... does not know difference between beginning and middle. (zeros interpreted here as vectors)
    if debug:  print("Subdiagonals (first array is lowest band, last array is diagonal band):\n",b)
    # execute fixed number of recursions to populate b and (possibly) v
    if block_action is None:
       for i in range(n):
            _iteration(H,v,b,reorthonormalize)         # In practice, there are too many cases where reorthonormalization is necessary for stability
            if debug and save:  print("Deviation from ON:", v[B:].orthonormality())
            if not save:  del v[0]                     # Release oldest vector, if not saving.
    else:
       bigN = n // block_action		# How many times to produce block_action number of new vectors
       liln = n - bigN*block_action	# The remainder to bring the total number of vectors up to n
       for i in range(bigN):
           _block_iteration(H,v,b,block_action,reorthonormalize)
           if debug and save:  print("Deviation from ON:", v[B:].orthonormality())
           if not save:  del v[:block_action]              # Release oldest vectors, if not saving.
       if liln>0:
           _block_iteration(H,v,b,liln,reorthonormalize)
           if debug and save:  print("Deviation from ON:", v[B:].orthonormality())
           if not save:  del v[:liln]                      # Release oldest vectors, if not saving.
    if save:  del v[:B]                                    # If not saved, all implicitly destroyed on return, but these are meaningless.
    if debug:  print("Vector addresses:\n",list(v))
    # populate a band diagonal matrix from b
    if autocomplete:     # Perhaps extend dimensionality by B and approximate lower-right block
        n += B
        for b_i in b[:B]:  b_i += zeros
        b[B] += b[B][-B:]
    pHp = numpy.matrix([[H.field(0)]*n]*n)
    for i in range(n):  pHp[i,i] = b[B][i]
    for i in range(1,B+1):
        for j in range(n-i):
            if debug:  print("pHp[{r},{c}] = pHp[{r},{c}]* = {b}".format(r=j+i,c=j,b=b[B-i][i+j]))
            pHp[j+i,j] = b[B-i][i+j]
            pHp[j,j+i] = conjugate(b[B-i][i+j])
    if debug:  print("Finished projection:\n",pHp)
    # return the band diagonal matrix, and, perhaps, the Lanczos vectors
    if save:  value = (pHp,v)
    else:     value =  pHp
    return value

def projection_eval_decomp(H, vecs, dim=10, block_action=None, autocomplete=False):
    """\
    Given a list (iterable type) of B orthonormal starting vectors in vecs, this function builds a
    square matrix projection of the linear operator H with dimension B*dim x B*dim, and it returns
    the eigenvalues of that projection along with all of the "short" eigenvectors and the "long"
    resoltions of the underlying Lanczos basis.

    Orthonormality of the input vecs is enforced by an in place Gram-Schmidt of the vectors therein.
    If block_action is set to an integer, it creates the Lanczos vectors by acting on block_action
    number of vectors at a time (or less in the case of a nonzero remainder to n/block_action ... 
    block_action cannot be greater than B, and if it is, it is quietly set equal to B at the beginning
    of the algorithm).

    If the autocomplete option is used, then a projection
    of dimension B larger than known explicitly is returned, where the bottom-right matrix block
    is approximated using the diagonals from the block that precedes it.  
    """
    debug = False	# Set to true to turn on verbose printing
    B = len(vecs)	# block size
    if debug:  print("Entering lanczos.projection_eval_decomp with B={}, dim={}, block_action={}, and autocomplete={}.".format(B,dim,block_action,autocomplete))
    # Lanczos projection of H with fixed dimension dim*B, and the vector set spanning the space of the projection
    # Use reorthonormalize=True because we explicitly want vectors and so need to be able to trust them.
    pHp,lanczos_vecs = projection(H, vecs, dim*B, block_action=block_action, save=True, autocomplete=autocomplete, reorthonormalize=True)	# "fixes" non-ON vecs or too large block_action
    # Diagonalize the projection and sort the eigenpairs
    eigen_vals,eigen_vecs = numpy.linalg.eigh(pHp)
    eigen_vals,eigen_vecs = zip(*sorted(zip( eigen_vals, eigen_vecs.T.tolist() ),key=lambda p: p[0]))
    return eigen_vals,eigen_vecs,lanczos_vecs

def projection_lowest_eigen(H, vecs, dim=10, block_action=None, autocomplete=False):
    """\
    Given a list (iterable type) of B orthonormal starting vectors in vecs, this function builds a
    square matrix projection of the linear operator H with dimension B*dim x B*dim, and it returns
    the lowest B eigenvalues of the projection and the "long" resolutions of the eigenvectors
    (so B approximate eigenvectors of H).

    Orthonormality of the input vecs is enforced by an in place Gram-Schmidt of the vectors therein.
    If block_action is set to an integer, it creates the Lanczos vectors by acting on block_action
    number of vectors at a time (or less in the case of a nonzero remainder to n/block_action ... 
    block_action cannot be greater than B, and if it is, it is quietly set equal to B at the beginning
    of the algorithm).

    If the autocomplete option is used, then a projection
    of dimension B larger than known explicitly is returned, where the bottom-right matrix block
    is approximated using the diagonals from the block that precedes it.  
    """
    debug = False	# Set to true to turn on verbose printing
    B = len(vecs)	# block size
    if debug:  print("Entering lanczos.projection_lowest_eigen with B={}, dim={}, block_action={}, and autocomplete={}.".format(B,dim,block_action,autocomplete))
    # Get the eigenvalue decomposition of the projection and the lanczos basis vectors
    eigen_vals, eigen_vecs, lanczos_vecs = projection_eval_decomp(H, vecs, dim, block_action=block_action, autocomplete=autocomplete)		# "fixes" non-ON vecs or too large block_action
    # Build the lowest B eigenvectors in the orginal full space
    augment = [] if autocomplete else [0]*B
    new_vecs = [ lanczos_vecs.deproject(eigen_vec+augment) for eigen_vec in eigen_vecs[:B] ]
    new_vals = eigen_vals[:B]
    return new_vals, new_vecs

def _lowest_eigen(H, vals, vecs, dim, block_action, autocomplete, converge_vectors, thresh):
    """\
    This is a helper function that extracts one (or more, read on) eigenpairs from the linear operator H,
    given B orthonormal starting vectors in vecs, which is a vector_set.  Orthonormality of the input vecs
    is enforced by an in place Gram-Schmidt of the vectors therein at the beginning of the algorithm.
    vals is an array of approximate eigenvalues for each of the vecs, assuming they approximate eigenvectors.
    If they do not, a list of B "infinities" should be given as input.  Differences relative to these
    values are used to judge convergence after the first iteration.

    Each iteration consists of building a Lanczos projection of dimension B*dim x B*dim
    of the operator whose eigenpairs are to be extracted (so, dim is the dimension of the Krylov space
    associated with each starting vector of each iteration).  If block_action is set to an integer,
    then the operator H is acted simultaneously on that number of vectors, block_action cannot be greater
    than B; if it is larger than B it is quietly reset to B.

    If the autocomplete option is used, then a projection of dimension B larger than known explicitly is
    diagonalized in each iteration, where the bottom-right matrix block is approximated using the diagonals
    from the block that precedes it.

    Convergence is met, if, after a given iteration, the change in one of the eigenvalues is smaller than
    thresh for ANY of the vectors, unless converge_vectors is True, in which case the largest projection of
    a vector outside the space spanned by the previous set should have a norm less than thresh for convergence.

    The eigenvalues, eigenvectors and respective error metrics for each vector are returned;
    at least one of the errors will be less than thresh.  These are not sorted, but each is in the same order.  
    vecs is returned as a vector_set.
    """
    debug = True			# Set to true to turn on verbose printing
    B = len(vecs)			# block size
    if debug:  print("Entering lanczos._lowest_eigen with B={}, vals={}, dim={}, block_action={}, autocomplete={}, and thresh={} ({}).".format(B,vals,dim,block_action,autocomplete,thresh,"vectors" if converge_vectors else "values"))
    # Repeatedly diagonalize projections until convergence
    errors = [float("inf")]*B
    t0 = time.time()
    iteration = 0
    while min(errors)>thresh:
        # Build the lowest B eigenvectors of the dim*B(+B?) projection of H
        new_vals, new_vecs = projection_lowest_eigen(H, vecs, dim, block_action=block_action, autocomplete=autocomplete)	# "fixes" non-ON vecs (drift with iterations) or too large block_action
        # Compute the convergence tests depending on whether we are testing vectors or values
        if converge_vectors:		# compute the projection of each outside the space spanned by the previous set
            P = 1 - vecs.projection_opr()
            errors = [(Vi|P|Vi) for Vi in new_vecs]
        else:				# compute the eigenvalue change
            errors = [abs(Enew-Eold) for Enew,Eold in zip(new_vals,vals)]
	# reset the guess
        vecs = vector_set(H.space,new_vecs)
        vals = list(new_vals)		# make list or it is a tuple that does not support deletion in call from lowest_eigen
        # admin
        t1 = time.time()
        if debug:
            print("Lanczos iteration {}:  Eigenvalues = {}\nerrors = {}".format(iteration,vals,errors))
            print('cycle time =', t1 - t0)
        t0 = t1
        iteration += 1
    # Give the resulting eigensolutions and their distances from convergence
    return vals,vecs,errors

def lowest_eigen(H, v, thresh, dim=10, num=None, block_action=None, autocomplete=False, converge_vectors=False):
    """\
    This function uses the iterative Lanczos method to extract the lowest num eigenpairs from a linear operator H,
    given B orthonormal initial guess vectors supplied in a list v (actually any iterable sequence).  Orthonormality
    of the input v is enforced by an in place Gram-Schmidt of the vectors therein at the beginning of the algorithm.

    The input vectors must have non-zero overlap with the targeted eigenvectors.
    num cannot be greater than B, but it can be smaller.  There are uses for having num greater than B, but such
    approaches should probably be encased in a different algorithm, the danger with restarted algorithms is that
    the periodic truncation refines the symmetries of the vectors (if there is any) and higher vectors with symmeties
    other than those in the block might be expunged.  It is better to have at least one input per output.

    If the autocomplete option is used, then a projection of dimension B larger than known explicitly is
    diagonalized in each iteration, where the bottom-right matrix block is approximated using the diagonals
    from the block that precedes it.

    Each iteration consists of building a Lanczos projection of dimension B*dim x B*dim
    of the operator whose eigenpairs are to be extracted (so, dim is the dimension of the Krylov space
    associated with each starting vector of each iteration).  If block_action is set to an integer,
    then the operator H is acted simultaneously on that number of vectors, which can be more efficient
    if there is a lot of logic involved in the operator action.  block_action cannot be greater than B;
    if it is larger than B it is quietly reset to B.

    At the moment when any given vector is converged, it is removed from the calculation and projected
    out of the space of consideration, in order to preserve stability.
    Overall convergence is met, if, after a given iteration, the change in each of the eigenvalues is smaller than
    thresh, unless converge_vectors is True, in which case the largest projection of
    the vectors outside the space spanned by the previous set should have a norm less than thresh for convergence.
    Converging the vectors is preferable
    when possible, since the criterion is independent of the scale of the matrix eigenvalues or their
    differences.  However, in cases where only a subset of a degenerate set are to be returned, then these
    are ill-defined and the algorithm never converges.  On the other hand, with the energy criterion,
    these ill-defined vectors are quietly returned, since this is a difficult condition to be sure of.
    """
    debug = True
    Bo = len(v)					# the initial block size (will shrink as vectors converge)
    if (num is None) or (num>Bo):  num = Bo	# if not specified, assume we want one eigenvector per starting vector
    if debug:  print("Entering lanczos.lowest_eigen with Bo={}, dim={}, num={}, block_action={}, autocomplete={}, and thresh={} ({}).".format(Bo,dim,num,block_action,autocomplete,thresh,"vectors" if converge_vectors else "values"))
    frozen_vals = []			# store converged eigenvalues
    frozen_vecs = vector_set(H.space)	# store converged eigenvectors
    vals        = [float("inf")]*Bo	# initial guess at eigenvalues (guaranteed to force second iteration if eigenvalues tested)
    vecs        = vector_set(H.space,v)	# promote to a vector set
    B = Bo
    while B>Bo-num or min(vals or [float("inf")])<max(frozen_vals or [-float("inf")]):				# vecs will constantly be replaced, when desired converged vectors are removed, we stop
        P = 1 - frozen_vecs.projection_opr()			# P projects out frozen vectors
        vecs = vector_set(H.space, [ P|Vi for Vi in vecs ])	# If we make sure they are projected out here ... (will be nicer when P|vecs works, also _lowest_eigen enforces ON)
        vals,vecs,errors = _lowest_eigen(P|H, vals, vecs, dim, block_action, autocomplete, converge_vectors, thresh) # ... then we only need P|H here and not P|H|P (b/c projectors are idempotnent, also block_action "fixed" at lower level)
        for n,error in reversed(list(enumerate(errors))):	# loop backwards so that deletions do not change indexes looped over later in time (enumerate is not itself reversible)
           if error<thresh:					# MOVE any converged vectors from vecs to frozen_vecs
                frozen_vals += [vals[n]]
                frozen_vecs += [vecs[n]]
                del vals[n]
                del vecs[n]
        B = len(vecs)
    return sorted(zip(frozen_vals,frozen_vecs), key=lambda p: p[0])	# return as list of eigenpair tuples sorted low to high
