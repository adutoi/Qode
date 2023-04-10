#    (C) Copyright 2016 Yuhong Liu
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
import sys

def ab_spin(num_e,multi):
    # number of electrons: num_e
    # multiplicity: multi
    # Return numbers of alpha and beta electrons.
        
    if ((num_e > 1) and (multi > 0)): 
        single_e = multi - 1
        if (num_e - single_e) < 0 or (num_e - single_e)%2 != 0:
            print("Invalid Charge/Electron Configuration")
            sys.exit(1)
        beta_e = (num_e - single_e) // 2
        alpha_e = beta_e + single_e
        
    elif (num_e == 1) and (multi == 2):
        alpha_e = 1
        beta_e = 0

    else:
        print("Invalid Charge/Electron Configuration")
        sys.exit(1)
        
    return alpha_e,beta_e
