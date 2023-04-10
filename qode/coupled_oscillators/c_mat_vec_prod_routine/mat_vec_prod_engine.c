/*   (C) Copyright 2018 Anthony D. Dutoi
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
    This diag_act_on_vec function takes a list of eigenvalues for ONE FRAGMENT, 
    and act those numbers on the original vector.
*/
int diag_act_on_vec(
                    int       total_ld        ,  /* Leading Dimension of the vector ( H Matrix as well ) */
                    double  *eigval_list     ,  /* ONE of the eigenvalue list (1D array) */
                    int       num_state       ,  /* number of states of current fragment  */
                    int       jumps           ,  /* Least jump if current fragment state-number plus one */
                    double  *orig_vec        ,  
                    double  *new_vec
                    )
{    
    int x = 0;
    while (x < total_ld)  /* They are looped many times to fill out the whole matrix dimension */
        {
            for (int state_i = 0 ; state_i < num_state ; state_i++ ) /* increase state number */
            {
                for (int jump_i = 0 ; jump_i < jumps; jump_i++)     /* Inside the SAME JUMP, Fill in the SAME VALUE */
                {
                    new_vec[x] += orig_vec[x] * eigval_list[state_i];
                    x++;
                }
            }
        }
    
    
    return 0;
}




/*
    coupling_act_on_vec function takes a pair of coupled Fragments and act this coupling matrix on vector orig_vec.
    It saves the changes to another vector ( new_vec ) which is requried at the time when this function got called.
    
*/
int coupling_act_on_vec(
                        int num_state1,
                        int num_state2,         /* shift in coupling mat = num_state2 */
                        int num1_jump,
                        int num2_jump,
                        int total_ld,           /* Leading Dimension of the vector ( H Matrix as well ) */
                        double *coupling_mat,   /* coupling matrix for this pair                        */
                        int left_limit,         /* loop limits */
                        int mid_limit,          /* loop limits */
                        int right_limit,        /* loop limits */
                        int left_jump,
                        int mid_jump,
                        int right_jump,
                        double *orig_vec,
                        double *new_vec
                        )
{
    int x, y;
    int coupling_mat_dim = num_state1 * num_state2;
    
    for (int c1_s1 = 0; c1_s1 < num_state1; c1_s1++)
    {
        for (int c2_s1 = 0; c2_s1 < num_state2; c2_s1++)
        {
            for (int c1_s2 = 0; c1_s2 < num_state1; c1_s2++)
            {
                for (int c2_s2 = 0; c2_s2 < num_state2; c2_s2++)
                {
                    double coupling_mat_value = coupling_mat[ (c1_s1*num_state2 + c2_s1) * coupling_mat_dim + c1_s2*num_state2 + c2_s2];
                    
                    for (int left_i = 0; left_i < left_limit; left_i++)
                    {
                        for (int mid_i = 0; mid_i < mid_limit; mid_i++)
                        {
                            for (int right_i = 0; right_i < right_limit; right_i++)
                            {
                                x = left_i * left_jump + \
                                    c1_s1  * num1_jump + \
                                    mid_i  * mid_jump  + \
                                    c2_s1  * num2_jump + \
                                    right_i;
                                
                                y = left_i * left_jump + \
                                    c1_s2  * num1_jump + \
                                    mid_i  * mid_jump  + \
                                    c2_s2  * num2_jump + \
                                    right_i;
                                
                                new_vec[x] += orig_vec[y] * coupling_mat_value;
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}



