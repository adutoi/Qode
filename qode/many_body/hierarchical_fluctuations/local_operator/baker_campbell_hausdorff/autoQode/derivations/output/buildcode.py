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
import os
from multiprocessing import Pool
from c_generator import add_obj_call_code
import pickle


def wrapper(inputfile):
	print(inputfile)
	prefix     = inputfile.split('.')[0]
	expression = pickle.load(open(inputfile,'rb'))
	add_obj_call_code(expression, prefix)


allp = [ item for item in os.listdir() if item.endswith('.p') ]

# allp  = [ 'Fv_Ex.p' ]
# allp  = [ 'FvFv_Ex.p' ]
# allp = ['FoFv_ExEx.p']
# allp = ['DxDx_ExEx.p']
# allp = ['Dx_ExEx.p']

# allp = ['ExFoFo_Ex.p']


for item in allp:
	wrapper(item)

# myPool = Pool(1)
# myPool.map(wrapper, allp)
