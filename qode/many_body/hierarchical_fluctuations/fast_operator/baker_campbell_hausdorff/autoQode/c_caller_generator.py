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
import sys
import os



def read_void_func(header_file):
	f = open('derivations/output/' + header_file, 'r')
	lines = f.readlines()
	f.close()
	void_funcs = [item for item in lines if "void" in item]
	return void_funcs


if __name__ == "__main__":
	all_headers = [ item for item in os.listdir('derivations/output') if item.endswith('.h') ]
	# print(all_headers)

	file_buf = "import ctypes as ct\nimport numpy as np\n"

	lib_path = os.getcwd() + "/derivations/output"

	for item in all_headers:
		print(item, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		file_main_name = item.split('.')[0] 
		wrapper_file   = file_main_name + '.py'
		print("WRAPPER NAME =", wrapper_file)
		output_buf  = "def %s( X_new, block_X, block_T, No, Nv , nthd ):\n" %(file_main_name)
		output_buf += "\tX_new_dict = X_new.__dict__\n\tc_handler = ct.cdll.LoadLibrary('%s/%s.so')\n" %(lib_path,file_main_name)
		# output_buf += "\tprint('CALLING %s.so')\n" %(file_main_name)
		output_buf += "\tc_double_ptr = ct.POINTER(ct.c_double)\n"
		output_buf += "\tX_ptr = block_X.ctypes.data_as(c_double_ptr)\n\tT_ptr = block_T.ctypes.data_as(c_double_ptr)\n"
		all_funcs  = read_void_func(item)
		for line in all_funcs:
			parts = line.split()
			# print(parts)
			for part in parts:
				if 'block' in part:
					key = part.split('_')[0]
			# print(parts, key)
			output_buf += '\t#\n\t%s_ptr = X_new_dict["%s"].ctypes.data_as(c_double_ptr)\n' %(key,key)
			output_buf += '\tc_handler.%sX_ptr, T_ptr, ct.c_int(No), ct.c_int(Nv), %s_ptr, ct.c_int(nthd) )\n' %(parts[1], key)
		output_buf += '\treturn X_new\n'
		# print(output_buf)
		file_buf += output_buf + '\n'

	# print(file_buf)
	f = open("c_caller.py",'w')
	f.write(file_buf)
	f.close()
