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

def read_wrapper_func(header_file):
	f = open('derivations/output/' + header_file, 'r')
	lines = f.readlines()
	f.close()
	void_funcs = [item for item in lines if "wrapper" in item]
	return void_funcs[0]



all_headers = [ item for item in os.listdir('derivations/output') if item.endswith('.h') ]
# all_headers = [ 'DxDx_ExEx.h' ]


# file_buf = "import ctypes as ct\nimport numpy as np\n"



# lib_path = os.environ['PYTHONPATH']
# if ':' in lib_path:
	# raise AssertionError

lib_path = os.getcwd() + "/derivations/output"

wrapper_func_name = "import ctypes as ct\nimport numpy as np\n\ndef compute_next_X( prev_X_obj, prev_T_obj, new_X_obj, rec_num_states, n_threads ):\n"
wrapper_func_body = "\tNo = ct.c_longlong(1)\n\tNv = np.zeros(len(rec_num_states), dtype=np.int64)\n\tfor i in range(len(rec_num_states)):\n\t\tNv[i] = rec_num_states[i] - 1\n"
# wrapper_func_body += "\tNv_c = (ct.c_longlong * len(Nv))()\n\tfor i in range(len(Nv)):\n\t\tNv_c[i] = Nv[i]\n\tp_Nv_c = ct.cast(Nv_c, ct.POINTER(ct.c_longlong))\n"
wrapper_func_body += "\tp_Nv_c = Nv.ctypes.data_as( ct.POINTER(ct.c_longlong))\n"
wrapper_func_body += "\tNMol   = ct.c_longlong(len(rec_num_states))\n\tnthd = ct.c_longlong(n_threads)\n"

for item in all_headers:
	prefix = item.split('.')[0]
	print("Running %s Generation" %(prefix))
	# print(prefix, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	wrapper_func_body += '\t#\n\t#\n\t#\n\t# %s Recipe\n\t#\n' %(prefix)
	wrapper_line = read_wrapper_func(item)
	X_block_name, T_block_name = prefix.split('_')
	# print("X =",X_block_name, "T =",T_block_name)
	wrapper_func_body += '\tif prev_X_obj.%s is not None and prev_T_obj.%s is not None:\n' %(X_block_name, T_block_name)
	# wrapper_func_body += '\t\tprint("[%s, %s] -----------------------------")\n' %(X_block_name,T_block_name)
	wrapper_func_body += '\t\tprev_X_dim = ct.c_longlong( prev_X_obj.%s_dim )\n\t\tprev_T_dim = ct.c_longlong( prev_T_obj.%s_dim )\n' %(X_block_name,T_block_name)
	# wrapper_func_body += '\t\tprint("prev_X_dim =", prev_X_obj.%s_dim , "prev_T_dim =", prev_T_obj.%s_dim)\n' %(X_block_name,T_block_name)
	wrapper_func_body += '\t\tblock_dim  = ct.c_longlong(0)\n'
	# wrapper_func_body += '\t\tprint("Running Recipe: %s")\n' %(prefix)
	# wrapper_func_body += '\tif True:\n'
	# wrapper_func_body += '\t\tprint("Running Recipe: %s")\n' %(prefix)
	type_stripped = [ item for item in wrapper_line.split() if "double" not in item and "int64_t" not in item and "void" not in item ]
	# print(type_stripped)
	# print(type_stripped[10:-5:])
	updating_blocks = type_stripped[10:-5]
	if len(updating_blocks) %3 != 0: # Must be a multiple of 3
		raise AssertionError

	wrapper_func_body += '\t\tprev_X_%s = prev_X_obj.%s.ctypes.data_as( ct.POINTER( ct.c_double ) )\n' %(X_block_name,X_block_name)
	wrapper_func_body += '\t\tX_starters = prev_X_obj.%s_starters.ctypes.data_as( ct.POINTER( ct.c_longlong ) )\n' %(X_block_name)
	wrapper_func_body += '\t\tX_dim = ct.c_longlong( prev_X_obj.%s_int_dim )\n' %(X_block_name)

	wrapper_func_body += '\t\tprev_T_%s = prev_T_obj.%s.ctypes.data_as( ct.POINTER( ct.c_double ) )\n' %(T_block_name,T_block_name)
	wrapper_func_body += '\t\tT_starters = prev_T_obj.%s_starters.ctypes.data_as( ct.POINTER( ct.c_longlong ) )\n' %(T_block_name)
	wrapper_func_body += '\t\tT_dim = ct.c_longlong( prev_T_obj.%s_int_dim )\n' %(T_block_name)

	
	wrapper_func_body += '\t\t%s_handle = ct.cdll.LoadLibrary("%s/derivations/output/%s.so")\n' %(prefix,os.getcwd(), prefix)

	# new_block_code = ""
	for i in range(0,len(updating_blocks),3):
		new_block_name = updating_blocks[i+1].split('_starter')[0]
		# new_block_code += ''
		wrapper_func_body += '\t\tnew_%s_block = new_X_obj.%s.ctypes.data_as( ct.POINTER( ct.c_double ) )\n' %(new_block_name, new_block_name)
		wrapper_func_body += '\t\t%s_starters  = new_X_obj.%s_starters.ctypes.data_as( ct.POINTER( ct.c_longlong ) )\n' %(new_block_name, new_block_name)
		wrapper_func_body += '\t\t%s_dim = ct.c_longlong( new_X_obj.%s_int_dim )\n' %(new_block_name, new_block_name)

	# print("***********+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

	calling_line = "\t\t%s_handle." %(prefix)
	for item in type_stripped[:-5]:
		calling_line += item + ' '
	calling_line += ' No, p_Nv_c, NMol, nthd)\n'
	wrapper_func_body += calling_line

file_buf = wrapper_func_name + wrapper_func_body
# print(file_buf)



f = open("c_caller.py",'w')
f.write(file_buf)
f.close()
