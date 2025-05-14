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
All Hamiltonian terms are defined in h01.py to h13.py.
First Run python3.3 gen_cisd_code_caller.py to derive and build the code for CISD. The new code is named normal_order_h_cisd_py_generated.py which is NOT CHECKED IN.
Then  Run python3.3 generated_code_caller.py to compute CISD energy for a saved F_mat.npy and V_mat.npy
