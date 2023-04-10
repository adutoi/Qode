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
from summation import Sum


# index_dict = { 'i':'No', 'j':'No', 'k':'No', 'l':'No', 'm':'No', 'n':'No', 'o':'No', \
# 	           'a':'Nv', 'b':'Nv', 'c':'Nv', 'd':'Nv', 'e':'Nv', 'f':'Nv', 'g':'Nv' }


def index_dict(index_with_mol):
	# index_with_mol is e.g. ['i', 'P']  ( 'i' is index; 'P' is molecule number )
	idx, mol = index_with_mol
	if idx in ['i','j','k','l','m','n','o','p']:
		return 'No'
	elif idx in ['a','b','c','d','e','f','g','h']:
		return 'Nv[%s]' %(mol)
	else:
		raise AssertionError

def form_index(index_with_mol):
	return "%s_%s" %(index_with_mol[0], index_with_mol[1])


def two_index_to_position( list_of_index ):
	if len(list_of_index) != 2:
		raise AssertionError
	#
	#
	piece2 = " + " + form_index(list_of_index[1])
	piece1 =         form_index(list_of_index[0]) + " * " + index_dict(list_of_index[1])
	return piece1 + piece2

def four_index_to_position( list_of_index ):
	if len(list_of_index) != 4:
		raise AssertionError

	#
	#  Assume this indexing scheme:
	#    Absolute Position in Array = idx0 * ( num_orb_1 * num_orb_2 * num_orb_3 ) + idx1 * ( num_orb_2 * num_orb_3 ) + idx2 * num_orb_3 + idx3 
	#  Generate this string...
	#
	piece4 = " + " + form_index(list_of_index[3])
	piece3 = " + " + form_index(list_of_index[2]) + " * " + index_dict(list_of_index[3])
	piece2 = " + " + form_index(list_of_index[1]) + " * " + index_dict(list_of_index[2]) + " * " + index_dict(list_of_index[3]) 
	piece1 =         form_index(list_of_index[0]) + " * " + index_dict(list_of_index[1]) + " * " + index_dict(list_of_index[2]) + " * " + index_dict(list_of_index[3])
	# print(piece1 + piece2 + piece3 + piece4)
	return piece1 + piece2 + piece3 + piece4



def six_index_to_position( list_of_index ):
	if len(list_of_index) != 6:
		raise AssertionError

	piece6 = " + " + form_index(list_of_index[5])
	piece5 = " + " + form_index(list_of_index[4]) + " * " + index_dict(list_of_index[5])
	piece4 = " + " + form_index(list_of_index[3]) + " * " + index_dict(list_of_index[4]) + " * " + index_dict(list_of_index[5]) 
	piece3 = " + " + form_index(list_of_index[2]) + " * " + index_dict(list_of_index[3]) + " * " + index_dict(list_of_index[4]) + " * " + index_dict(list_of_index[5])
	piece2 = " + " + form_index(list_of_index[1]) + " * " + index_dict(list_of_index[2]) + " * " + index_dict(list_of_index[3]) + " * " + index_dict(list_of_index[4]) + " * " + index_dict(list_of_index[5])
	piece1 =         form_index(list_of_index[0]) + " * " + index_dict(list_of_index[1]) + " * " + index_dict(list_of_index[2]) + " * " + index_dict(list_of_index[3]) + " * " + index_dict(list_of_index[4]) + " * " + index_dict(list_of_index[5])

	return piece1 + piece2 + piece3 + piece4 + piece5 + piece6



def index_to_position( list_of_index ):
	if len(list_of_index) == 4:
		return four_index_to_position(list_of_index)
	elif len(list_of_index) == 2:
		return two_index_to_position(list_of_index)
	elif len(list_of_index) == 6:
		return six_index_to_position(list_of_index)
	elif len(list_of_index) == 0:
		return ''
	else:
		raise AssertionError


def append_to_file(file_name, output):
	f = open(file_name, 'a')
	f.write(output)
	f.close()




def one_term_gen(outmost_char_bank, index_sets, summand, func_num, prefix):
	print("ONE-TERM CASE IS NOT IMPLEMENTED FOR NOW!")
	raise AssertionError




def loop_condition_analyzer(sum_indices, loop_upper, loop_lower):
	# print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	# print("INDICES =", sum_indices)
	# print("UPPERS  =", loop_upper)
	# print("LOWERS  =", loop_lower)
	occ_idx = ['i','j','k','l','m','n','o','p']
	vrt_idx = ['a','b','c','d','e','f','g','h']
	num_loops = len(sum_indices)
	occ_loops = []
	vrt_loops = []
	for i in range(num_loops):
		if sum_indices[i][0] in occ_idx:
			# print(sum_indices[i], loop_lower[i], loop_upper[i])
			occ_loops += [[ sum_indices[i], loop_lower[i], loop_upper[i] ]]
		elif sum_indices[i][0] in vrt_idx:
			# print(sum_indices[i], loop_lower[i], loop_upper[i])
			vrt_loops += [[ sum_indices[i], loop_lower[i], loop_upper[i] ]]
		else:
			raise AssertionError
	
	continue_checking = True
	function_will_run = True 
	

	if occ_loops[0][1] == '0' and occ_loops[0][2] == 'No':
		# checking here... 
		for i in range(1, len(occ_loops)):
			if occ_loops[0][0] in occ_loops[i][1] or occ_loops[0][0] in occ_loops[i][2]:
				# print(occ_loops[0], 'in', occ_loops[i])
				continue_checking = False
				function_will_run = False
	else:
		raise AssertionError

	# if not dead ( continue_checking == True )
	i = 1
	while i < len(occ_loops) and continue_checking:
		if occ_loops[i][1] == '0' and occ_loops[i][2] == 'No':
			for j in range(i+1, len(occ_loops)):
				if occ_loops[i][0] in occ_loops[j][1] or occ_loops[i][0] in occ_loops[j][2]:
					# print(occ_loops[i], 'in', occ_loops[j])
					continue_checking = False
					function_will_run = False
		else:
			raise AssertionError
		i += 1
		

	# print('Will Run =', function_will_run)
	# print('----------------------------------------------------------')

	return function_will_run


def mol_indices_analyzer(operators):
	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	for item in operators:
		print(item)
	print('----------------------------------------------------------')
	return True





KEEP_TRIVIAL_FUNC = True
# KEEP_TRIVIAL_FUNC = False

# PRINT_VALUE_IN_C = True
PRINT_VALUE_IN_C = False

# FULL_CASE_PRINT = True
FULL_CASE_PRINT = False

def full_case_gen(outmost_char_bank, index_sets, summand, func_num, prefix):
	from base import copy

	if FULL_CASE_PRINT:
		print("FULL_CASE:")
		print("PREFIX =", prefix)
		print("OUTER CHAR BANK =", outmost_char_bank)



	X_suffix, T_suffix = prefix.split('_')
	#
	# OUTER SUMMATION PART:
	#
	#


	outer_sum_indices = []
	outer_loop_upper  = []
	outer_loop_lower  = []
	for outer_chars in outmost_char_bank:
		outer_sum_indices += [ outer_chars[0] ]
		outer_loop_upper  += [ outer_chars[1] ]
		outer_loop_lower  += [ '0' ]
		for relate, other in outer_chars[2:]:
			if relate == '<':
				outer_loop_upper[-1] = other
			elif relate == '>':
				outer_loop_lower[-1] = other + '+1'
			else:
				raise AssertionError
	if FULL_CASE_PRINT:
		print("OUTER LEVEL SUM OBJ ========================")
		print("outer_INDICES =", outer_sum_indices)
		print("outer_LOWERS  =", outer_loop_lower)
		print("outer_UPPERS  =", outer_loop_upper)
		print("=============================================")

	outer_func_name = "void compute_contraction_" + str(func_num) + '('
	outer_func_body = '{\nomp_set_dynamic(0);\nomp_set_num_threads(nthd);\n#pragma omp parallel for\n'

	for i in range( len(outer_sum_indices) ):
		outer_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(outer_sum_indices[i], outer_loop_lower[i], outer_sum_indices[i], outer_loop_upper[i], outer_sum_indices[i])



	#
	# NORMAL SUMMATION PART:
	#
	#
	op_indices = []
	for factor in summand.factors[1].factors:  
		op_indices += [factor.bottom, factor.top]
	char_bank = Sum._index_ordering(op_indices, index_sets, should_match=True)


	if FULL_CASE_PRINT:
		print("CHAR BANK",char_bank)


	for item in char_bank:
		if item[1] == 'Nv':
			item[1] = "%s[%s]" %(item[1] , item[0][-1])



	sum_indices = []
	loop_upper  = []
	loop_lower  = []
	for chars in char_bank:
		sum_indices += [ chars[0] ]
		loop_upper  += [ chars[1] ]
		loop_lower  += [ '0' ]
		for relate, other in chars[2:]:
			if relate == '<':
				loop_upper[-1] = other
			elif relate == '>':
				loop_lower[-1] = other + '+1'
			else:
				raise AssertionError

	if FULL_CASE_PRINT:
		print("NORMAL LEVEL SUM OBJ ========================")
		print("INDICES =", sum_indices)
		print("LOWERS  =", loop_lower)
		print("UPPERS  =", loop_upper)
		print("=============================================")


	for i in range( len(sum_indices) ):
		outer_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(sum_indices[i], loop_lower[i], sum_indices[i], loop_upper[i], sum_indices[i])
		




	target_block  = summand.factors[1].code()
	num_operators = len(target_block)

	# mol_indices_analyzer(target_block)
	# print(target_block)
	if num_operators == 3:
		outer_code = "%s*NMol*NMol + %s*NMol + %s" %(target_block[0][1][0][1], target_block[1][1][0][1], target_block[2][1][0][1])
	elif num_operators == 2:
		outer_code = "%s*NMol + %s" %(target_block[0][1][0][1], target_block[1][1][0][1])
	elif num_operators == 1:
		outer_code = "%s" %(target_block[0][1][0][1])
	else:
		raise AssertionError

	target_indices = []
	for i in range(num_operators):
		target_indices += target_block[i][1]

	block_code = ""
	for item in target_block:
		block_code += item[0]


	if FULL_CASE_PRINT:
		print("OUTER_CODE =", outer_code)
		print("TARGET_INDICES =",target_indices)

	outer_func_body += "#pragma omp atomic\n"
	outer_func_body += "%s_block[%s][%s] += " %(block_code, outer_code, index_to_position(target_indices))



	# print("Target Block =", target_block)
	#
	# INNER SUMMATION PART:
	#
	#
	inner_sum = summand.factors[0]
	inner_index_sets, inner_summand = inner_sum._pass_through_helper(copy, [])
	
	X = inner_summand.factors[-2]
	T = inner_summand.factors[-1]

	inner_op_indices = []
	possible_minus = inner_summand.factors[0] 	# possibly the same as X (ok, b/c just testing that it is Coeff)
	
	inner_coeff_indices = []					# all the index pairs of the Coeffs, with lowest order indices at the end
	for pair in range(max(len(X.top),len(T.top))):		# bottom and top of each have same length
		p = -1-pair
		if pair<len(T.top):     inner_coeff_indices = [T.top[p]]    + inner_coeff_indices
		if pair<len(X.top):     inner_coeff_indices = [X.top[p]]    + inner_coeff_indices
		if pair<len(T.bottom):  inner_coeff_indices = [T.bottom[p]] + inner_coeff_indices
		if pair<len(X.bottom):  inner_coeff_indices = [X.bottom[p]] + inner_coeff_indices
	inner_char_bank = Sum._index_ordering(inner_coeff_indices, inner_index_sets)

	inner_terms = inner_summand
	# print("TERMS =", end='  ')
	# for inner_term in inner_terms.factors:
	# 	term_code = inner_term.code()
	# 	print(term_code, end='  ')
	# print(" ")

	inner_sum_indices = []
	inner_loop_upper  = []
	inner_loop_lower  = []

	if FULL_CASE_PRINT:
		print("INNER CHAR BANK",inner_char_bank)

	for item in inner_char_bank:
		if item[1] == 'Nv':
			item[1] = "%s[%s]" %(item[1] , item[0][-1])

	for chars in inner_char_bank:
		inner_sum_indices += [ chars[0] ]
		inner_loop_upper  += [ chars[1] ]
		inner_loop_lower  += [ '0' ]
		for relate, other in chars[2:]:
			if relate == '<':
				inner_loop_upper[-1] = other
			elif relate == '>':
				inner_loop_lower[-1] = other + '+1'
			else:
				raise AssertionError


	if FULL_CASE_PRINT:
		print("INNER SUM OBJ ========================")
		print("INDICES =", inner_sum_indices)
		print("LOWERS  =", inner_loop_lower)
		print("UPPERS  =", inner_loop_upper)
		print("=============================================")


	inner_func_name = "double compute_inner_contraction_" + str(func_num) + '('
	inner_func_body = "{\ndouble contraction_value = 0.0;\n"

	for i in range( len(inner_sum_indices) ):
		inner_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(inner_sum_indices[i], inner_loop_lower[i], inner_sum_indices[i], inner_loop_upper[i], inner_sum_indices[i])
	
	inner_func_body += "contraction_value += "



	print_buf = 'printf("value = '
	var_buf   = ''

	for inner_term in inner_terms.factors:
		term_code = inner_term.code()
		# print(term_code)
		if len(term_code[1]) > 0:
			if term_code[0] == "X":
				suffix = X_suffix
			elif term_code[0] == "T":
				suffix = T_suffix
			else:
				raise AssertionError
			if len(term_code[1]) == 2:
				inner_func_body += "%s_%s[%s][%s]*" %(term_code[0], suffix, term_code[1][0][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s][%s]" %(term_code[0], suffix, term_code[1][0][1], index_to_position(term_code[1]) )

			elif len(term_code[1]) == 4:
				inner_func_body += "%s_%s[%s*NMol + %s][%s]*" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s*NMol + %s][%s]" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], index_to_position(term_code[1]) )

			elif len(term_code[1]) == 6:
				inner_func_body += "%s_%s[%s*NMol*NMol + %s*NMol + %s][%s]*" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], term_code[1][4][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s*NMol*NMol + %s*NMol + %s][%s]" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], term_code[1][4][1], index_to_position(term_code[1]) )

			else:
				raise AssertionError

			if term_code[0] == "X":
				inner_func_name += " double** %s_%s," %(term_code[0], X_suffix)
				outer_func_name += " double** %s_%s," %(term_code[0], X_suffix)
			elif term_code[0] == "T":
				inner_func_name += " double** %s_%s," %(term_code[0], T_suffix)
				outer_func_name += " double** %s_%s," %(term_code[0], T_suffix)

		else:
			inner_func_body += term_code[0] + '*'
			print_buf += term_code[0] + 'x'
	print_buf = print_buf[:-1] +  '\\n"' + var_buf + ');\n'
	# print(print_buf)


	for inner_term in inner_terms.factors:
		term_code = inner_term.code()
		if len(term_code[1]) > 0:
			for item in term_code[1]:
				new_idx = "%s_%s" %(item[0], item[1])
				if new_idx in sum_indices and new_idx not in inner_func_name:
					inner_func_name += " int64_t %s_%s," %(item[0], item[1])
	for imol in outer_sum_indices:
		inner_func_name += " int64_t %s," %(imol)
	inner_func_name += ' int64_t No, int64_t* Nv, int64_t NMol )'

	# for inner_index in inner_sum_indices:
	# 	inner_func_name += " int %s," %(inner_index)


	inner_func_body = inner_func_body[:-1] + ';\n'

	if PRINT_VALUE_IN_C:
		inner_func_body += print_buf

	inner_func_body += '}\n' * (inner_func_body.count('{') -1 )
	# inner_func_body += 'if (contraction_value == 0.0)\n{\nprintf("Return Value = %lf\\n", contraction_value);\n}\n'
	inner_func_body += 'return contraction_value;\n}\n'


	type_stripped = [ item for item in inner_func_name.split() if 'int64_t' not in item and 'double' not in item ]
	# print(type_stripped)

	for item in type_stripped:
		outer_func_body += item + ' '

	outer_func_body += ";\n"
	outer_func_body += '}\n' * outer_func_body.count('{') 

	outer_func_name +=  " double** %s_block, int64_t No, int64_t* Nv, int64_t NMol, int64_t nthd )" %(block_code)

	if FULL_CASE_PRINT:
		print("------------------------------------------------------------")
		print(outer_func_name)
		print(outer_func_body)
		print("------------------------------------------------------------")
		print(inner_func_name)
		print(inner_func_body)
		print("------------------------------------------------------------")


	code_parts   = '\n' + inner_func_name + '\n'  + inner_func_body + '\n' + outer_func_name + '\n' + outer_func_body
	header_parts = inner_func_name + ';\n' + outer_func_name + ';\n'
	
	# print(code_parts)
	# print(header_parts)

	if KEEP_TRIVIAL_FUNC:
		append_to_file(prefix+'.c', code_parts)
		append_to_file(prefix+'.h', header_parts)
	else:
		if loop_condition_analyzer(sum_indices + inner_sum_indices, loop_upper + inner_loop_upper, loop_lower + inner_loop_lower):
           # and mol_indices_analyzer(sum_indices + inner_sum_indices, loop_upper + inner_loop_upper, loop_lower + inner_loop_lower):
			append_to_file(prefix+'.c', code_parts)
			append_to_file(prefix+'.h', header_parts)

























# TWO_MUL_PRINT = True
TWO_MUL_PRINT = False

def two_mult_gen(outmost_char_bank, index_sets, summand, func_num, prefix):
	X_suffix, T_suffix = prefix.split('_')

	if TWO_MUL_PRINT:
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print("TWO_MUL_CASE:")
		print("OUTER CHAR BANK =", outmost_char_bank)

	outer_sum_indices = []
	outer_loop_upper  = []
	outer_loop_lower  = []
	for outer_chars in outmost_char_bank:
		outer_sum_indices += [ outer_chars[0] ]
		outer_loop_upper  += [ outer_chars[1] ]
		outer_loop_lower  += [ '0' ]
		for relate, other in outer_chars[2:]:
			if relate == '<':
				outer_loop_upper[-1] = other
			elif relate == '>':
				outer_loop_lower[-1] = other + '+1'
			else:
				raise AssertionError
	if TWO_MUL_PRINT:
		print("OUTER LEVEL SUM OBJ ========================")
		print("outer_INDICES =", outer_sum_indices)
		print("outer_LOWERS  =", outer_loop_lower)
		print("outer_UPPERS  =", outer_loop_upper)
		print("=============================================")

	outer_func_name = "void compute_contraction_" + str(func_num) + '('
	outer_func_body = '{\nomp_set_dynamic(0);\nomp_set_num_threads(nthd);\n#pragma omp parallel for\n'

	for i in range( len(outer_sum_indices) ):
		outer_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(outer_sum_indices[i], outer_loop_lower[i], outer_sum_indices[i], outer_loop_upper[i], outer_sum_indices[i])

	op_indices = []
	for factor in summand.factors[1].factors:  
		op_indices += [factor.bottom, factor.top]
	char_bank = Sum._index_ordering(op_indices, index_sets, should_match=True)

	if TWO_MUL_PRINT:
		print("CHAR BANK =", char_bank)

	for item in char_bank:
		if item[1] == 'Nv':
			item[1] = "%s[%s]" %(item[1] , item[0][-1])

	sum_indices = []
	loop_upper  = []
	loop_lower  = []
	for chars in char_bank:
		sum_indices += [ chars[0] ]
		loop_upper  += [ chars[1] ]
		loop_lower  += [ '0' ]
		for relate, other in chars[2:]:
			if relate == '<':
				loop_upper[-1] = other
			elif relate == '>':
				loop_lower[-1] = other + '+1'
			else:
				raise AssertionError

	if TWO_MUL_PRINT:
		print("NORMAL LEVEL SUM OBJ ========================")
		print("INDICES =", sum_indices)
		print("LOWERS  =", loop_lower)
		print("UPPERS  =", loop_upper)
		print("=============================================")


	for i in range( len(sum_indices) ):
		outer_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(sum_indices[i], loop_lower[i], sum_indices[i], loop_upper[i], sum_indices[i])
		

	mult_obj1, mult_obj2 = summand.factors
	target_block = mult_obj2.code()

	if TWO_MUL_PRINT:
		print("Target Block =", target_block)

	num_operators = len(target_block)
	# print(target_block)
	if num_operators == 3:
		outer_code = "%s*NMol*NMol + %s*NMol + %s" %(target_block[0][1][0][1], target_block[1][1][0][1], target_block[2][1][0][1])
	elif num_operators == 2:
		outer_code = "%s*NMol + %s" %(target_block[0][1][0][1], target_block[1][1][0][1])
	elif num_operators == 1:
		outer_code = "%s" %(target_block[0][1][0][1])
	else:
		raise AssertionError

	target_indices = []
	for i in range(num_operators):
		target_indices += target_block[i][1]

	block_code = ""
	for item in target_block:
		block_code += item[0]

	# print(outer_code, target_indices, block_code)


	inner_terms = mult_obj1.factors
	if TWO_MUL_PRINT:
		print("TERMS =", end='  ')
		for inner_term in inner_terms:
			term_code = inner_term.code()
			print(term_code, end='  ')
		print(" ")

	outer_func_body += "#pragma omp atomic\n"
	outer_func_body += "%s_block[%s][%s] += " %(block_code, outer_code, index_to_position(target_indices))


	print_buf = 'printf("value = '
	var_buf   = ''

	for inner_term in inner_terms:
		term_code = inner_term.code()
		# print(term_code)
		if len(term_code[1]) > 0:
			if term_code[0] == "X":
				suffix = X_suffix
			elif term_code[0] == "T":
				suffix = T_suffix
			else:
				raise AssertionError
			if len(term_code[1]) == 2:
				outer_func_body += "%s_%s[%s][%s]*" %(term_code[0], suffix, term_code[1][0][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s][%s]" %( term_code[0], suffix, term_code[1][0][1], index_to_position(term_code[1]) )

			elif len(term_code[1]) == 4:
				outer_func_body += "%s_%s[%s*NMol + %s][%s]*" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s*NMol + %s][%s]" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], index_to_position(term_code[1]) )

			elif len(term_code[1]) == 6:
				outer_func_body += "%s_%s[%s*NMol*NMol + %s*NMol + %s][%s]*" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], term_code[1][3][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s*NMol*NMol + %s*NMol + %s][%s]" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], term_code[1][3][1], index_to_position(term_code[1]) )

			else:
				raise AssertionError
			if term_code[0] == "X":
				outer_func_name += " double** %s_%s," %(term_code[0], X_suffix)
			elif term_code[0] == "T":
				outer_func_name += " double** %s_%s," %(term_code[0], T_suffix)

		else:
			outer_func_body += term_code[0] + '*'
			print_buf += term_code[0] + 'x'

	print_buf = print_buf[:-1] +  '\\n"' + var_buf + ');\n'

	outer_func_body = outer_func_body[:-1] + ';\n'

	if PRINT_VALUE_IN_C:
		outer_func_body += print_buf

	outer_func_body += '}\n' * outer_func_body.count('{')
	outer_func_name +=  " double** %s_block, int64_t No, int64_t* Nv, int64_t NMol, int64_t nthd )" %(block_code)





	if TWO_MUL_PRINT:
		print("------------------------------------------------------------")
		print(outer_func_name)
		print(outer_func_body)
		print("------------------------------------------------------------")

	code_parts   = '\n' + outer_func_name + '\n' + outer_func_body
	header_parts = outer_func_name + ';\n'
	
	# print(code_parts)
	# print(header_parts)

	if KEEP_TRIVIAL_FUNC:
		append_to_file(prefix+'.c', code_parts)
		append_to_file(prefix+'.h', header_parts)
	else:
		if loop_condition_analyzer(sum_indices, loop_upper, loop_lower ):
			append_to_file(prefix+'.c', code_parts)
			append_to_file(prefix+'.h', header_parts)


















# THREE_COEFF_PRINT = True
THREE_COEFF_PRINT = False

def three_coeff_gen(outmost_char_bank, index_sets, summand, func_num, prefix):
	X_suffix, T_suffix = prefix.split('_')

	if THREE_COEFF_PRINT:
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print("3_COEFF_CASE:")
		print("OUTER CHAR BANK =", outmost_char_bank)

	outer_sum_indices = []
	outer_loop_upper  = []
	outer_loop_lower  = []
	for outer_chars in outmost_char_bank:
		outer_sum_indices += [ outer_chars[0] ]
		outer_loop_upper  += [ outer_chars[1] ]
		outer_loop_lower  += [ '0' ]
		for relate, other in outer_chars[2:]:
			if relate == '<':
				outer_loop_upper[-1] = other
			elif relate == '>':
				outer_loop_lower[-1] = other + '+1'
			else:
				raise AssertionError
	if THREE_COEFF_PRINT:
		print("OUTER LEVEL SUM OBJ ========================")
		print("outer_INDICES =", outer_sum_indices)
		print("outer_LOWERS  =", outer_loop_lower)
		print("outer_UPPERS  =", outer_loop_upper)
		print("=============================================")

	outer_func_name = "void compute_contraction_" + str(func_num) + '('
	outer_func_body = '{\nomp_set_dynamic(0);\nomp_set_num_threads(nthd);\n#pragma omp parallel for\n'

	for i in range( len(outer_sum_indices) ):
		outer_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(outer_sum_indices[i], outer_loop_lower[i], outer_sum_indices[i], outer_loop_upper[i], outer_sum_indices[i])



	X,T = summand.factors[-2], summand.factors[-1]
	coeff_indices = []					# all the index pairs of the Coeffs, with lowest order indices at the end
	for pair in range(max(len(X.top),len(T.top))):		# bottom and top of each have same length
		p = -1-pair
		if pair<len(T.top):     coeff_indices = [T.top[p]]    + coeff_indices
		if pair<len(X.top):     coeff_indices = [X.top[p]]    + coeff_indices
		if pair<len(T.bottom):  coeff_indices = [T.bottom[p]] + coeff_indices
		if pair<len(X.bottom):  coeff_indices = [X.bottom[p]] + coeff_indices
	char_bank = Sum._index_ordering(coeff_indices, index_sets)

	inner_terms = summand.factors

	if THREE_COEFF_PRINT:
		print("CHAR BANK",char_bank)
		print("Target Block = K_block[0]")
		print("TERMS =", end='  ')
		for inner_term in inner_terms:
			term_code = inner_term.code()
			print(term_code, end='  ')
		print(" ")

	for item in char_bank:
		if item[1] == 'Nv':
			item[1] = "%s[%s]" %(item[1] , item[0][-1])

	sum_indices = []
	loop_upper  = []
	loop_lower  = []
	for chars in char_bank:
		sum_indices += [ chars[0] ]
		loop_upper  += [ chars[1] ]
		loop_lower  += [ '0' ]
		for relate, other in chars[2:]:
			if relate == '<':
				loop_upper[-1] = other
			elif relate == '>':
				loop_lower[-1] = other + '+1'
			else:
				raise AssertionError



	if THREE_COEFF_PRINT:
		print("NORMAL LEVEL SUM OBJ ========================")
		print("INDICES =", sum_indices)
		print("LOWERS  =", loop_lower)
		print("UPPERS  =", loop_upper)
		print("=============================================")

	for i in range( len(sum_indices) ):
		outer_func_body += "for ( int64_t %s = %s; %s < %s; %s++ )\n{\n" %(sum_indices[i], loop_lower[i], sum_indices[i], loop_upper[i], sum_indices[i])
		
	outer_func_body += "#pragma omp atomic\n"
	outer_func_body += 'K_block[0][0] += '


	print_buf = 'printf("value = '
	var_buf   = ''

	inner_terms = summand.factors
	for inner_term in inner_terms:
		term_code = inner_term.code()
		# print(term_code)
		if len(term_code[1]) > 0:
			if term_code[0] == "X":
				suffix = X_suffix
			elif term_code[0] == "T":
				suffix = T_suffix
			else:
				raise AssertionError
			if len(term_code[1]) == 2:
				outer_func_body += "%s_%s[%s][%s]*" %(term_code[0], suffix, term_code[1][0][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s][%s]" %(term_code[0], suffix, term_code[1][0][1], index_to_position(term_code[1]) )
			elif len(term_code[1]) == 4:
				outer_func_body += "%s_%s[%s*NMol + %s][%s]*" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s*NMol + %s][%s]" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], index_to_position(term_code[1]) )
			elif len(term_code[1]) == 6:
				outer_func_body += "%s_%s[%s*NMol*NMol + %s*NMol + %s][%s]*" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], term_code[1][4][1], index_to_position(term_code[1]) )
				print_buf += " %lf x"
				var_buf += ",%s_%s[%s*NMol*NMol + %s*NMol + %s][%s]" %(term_code[0], suffix, term_code[1][0][1], term_code[1][2][1], term_code[1][4][1], index_to_position(term_code[1]) )
			else:
				raise AssertionError
			if term_code[0] == "X":
				outer_func_name += " double** %s_%s," %(term_code[0], X_suffix)
			elif term_code[0] == "T":
				outer_func_name += " double** %s_%s," %(term_code[0], T_suffix)

		else:
			outer_func_body += term_code[0] + '*'
			print_buf += term_code[0] + 'x'

	outer_func_body = outer_func_body[:-1] + ';\n'

	print_buf = print_buf[:-1] +  '\\n"' + var_buf + ');\n'

	if PRINT_VALUE_IN_C:
		outer_func_body += print_buf
	
	outer_func_body += '}\n' * outer_func_body.count('{')	
	outer_func_name +=  " double** K_block, int64_t No, int64_t* Nv, int64_t NMol, int64_t nthd )"


	# print(print_buf)

	if THREE_COEFF_PRINT:
		print("------------------------------------------------------------")
		print(outer_func_name)
		print(outer_func_body)
		print("------------------------------------------------------------")

	code_parts   = '\n' + outer_func_name + '\n' + outer_func_body
	header_parts = outer_func_name + ';\n'
	# print(code_parts)
	# print(header_parts)

	if KEEP_TRIVIAL_FUNC:
		append_to_file(prefix+'.c', code_parts)
		append_to_file(prefix+'.h', header_parts)	
	else:
		if loop_condition_analyzer(sum_indices, loop_upper, loop_lower ):
			append_to_file(prefix+'.c', code_parts)
			append_to_file(prefix+'.h', header_parts)





def module_wrapper(prefix):
	f = open(prefix + '.h' , 'r')
	contents = f.read()
	f.close()
	pieces = contents.split()
	all_blocks = [item for item in pieces if "_block" in item ]
	redundant_blocks = [item.split("_block")[0] for item in all_blocks]
	target_blocks = []
	for item in redundant_blocks:
		if item not in target_blocks:
			target_blocks += [item]
	input_blocks = prefix.split('_')
	# print("INPUT  BLOCKS =", input_blocks)
	# print("TARGET BLOCKS =", target_blocks)

	wrapper_func_name = "\nvoid %s_wrapper_func( int64_t prev_X_dim, int64_t prev_T_dim, int64_t block_dim, double* prev_X_%s, int64_t* X_starters, int64_t X_dim, double* prev_T_%s, int64_t* T_starters, int64_t T_dim, " %(prefix, input_blocks[0], input_blocks[1] )
	wrapper_func_body =  "{\ndouble** X_%s = assign_ptrs( prev_X_%s, X_starters, X_dim);\ndouble** T_%s = assign_ptrs( prev_T_%s, T_starters, T_dim);\n" %(input_blocks[0], input_blocks[0],input_blocks[1], input_blocks[1])
	# wrapper_func_body += "for (int i=0; i< prev_X_dim; i++)\n{\nprintf(\"prev_X_%s[%%d] = %%lf\\n\",i, prev_X_%s[i]);\n}\n" %(input_blocks[0],input_blocks[0])
	# wrapper_func_body += "for (int i=0; i< prev_T_dim; i++)\n{\nprintf(\"prev_T_%s[%%d] = %%lf\\n\",i, prev_T_%s[i]);\n}\n" %(input_blocks[1],input_blocks[1])

	
	for item in target_blocks:
		wrapper_func_name += "double* new_%s_block, int64_t* %s_starters, int64_t %s_dim, " %(item, item, item)
		wrapper_func_body += "double** %s_block = assign_ptrs( new_%s_block, %s_starters, %s_dim);\n" %(item,item,item,item)
	
	wrapper_func_name += " int64_t No, int64_t* Nv, int64_t NMol, int64_t nthd )"

	
	lines = contents.split('\n')
	for line in lines:
		if 'void' in line:
			type_stripped = line.split()
			newline = ""
			for item in type_stripped:
				if "double" not in item and "int64_t" not in item and "void" not in item:
					newline += item + ' '
			wrapper_func_body += newline + '\n'
			# wrapper_func_body += 'puts("%s");\n' %(type_stripped[1])
	
	wrapper_func_body += "}\n"

	# print("------------------------------------------------------------")
	# print(wrapper_func_name)
	# print(wrapper_func_body)
	# print("------------------------------------------------------------")
	append_to_file(prefix+'.h', wrapper_func_name + ';\n')
	append_to_file(prefix+'.c', wrapper_func_name + '\n' + wrapper_func_body)





def build_shared_lib(file_main_name):
	#
	# USE THIS WHEN ALL CODE GENERATIONS ARE DONE!
	#
	import os
	compile_cmd = "gcc -std=c99 -fopenmp -c -fPIC -O2 " + file_main_name + ".c -I. -o " + file_main_name + ".o"
	linking_cmd = "gcc -shared -fopenmp -O2 " + file_main_name  + ".o -o " + file_main_name + ".so"

	# print('_' * 80)
	print(compile_cmd)
	if os.system( compile_cmd ) != 0:
		raise RuntimeError

	if os.system( linking_cmd ) != 0:
		raise RuntimeError
	print(linking_cmd)
	print('_' * 80)


def add_obj_call_code(add_obj, prefix):
	import os
	if os.system("rm -f %s.c %s.h %s.o %s.so" %(prefix, prefix, prefix, prefix)) != 0:
		raise RuntimeError

	append_to_file(prefix+'.c', '#include <stdio.h>\n#include <stdlib.h>\n#include <stdint.h>\n#include "%s.h"\n#include "omp.h"\n' %(prefix))
	append_to_file(prefix+'.h', "\n")
	macros = \
'''


double** assign_ptrs(double* storage, int64_t* starters, int64_t dimension)
{
	// Make an array of double*, pointing to storage[ starters[i] ] locations
	double** p_arrays = (double**) malloc(sizeof(double*)*dimension);
	for (int64_t i = 0; i < dimension; i++)
	{
		p_arrays[i] = &storage[ starters[i] ];
	}
	return p_arrays;
}

'''

	append_to_file(prefix+'.c', macros)
	file_main_name = prefix
	# print("-----------------------------------------------------------------------------------------")
	# print("gcc -std=c99 -fopenmp -c -fPIC " + file_main_name + ".c -I. -o " + file_main_name + ".o")
	# print("gcc -shared -fopenmp " + file_main_name  + ".o -o " + file_main_name + ".so")
	# print("-----------------------------------------------------------------------------------------")
	ct = 0
	for item in add_obj.terms:
		# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		# print(item)
		item.code(ct, prefix)
		ct += 1
	module_wrapper(prefix)
	build_shared_lib(prefix)

'''
int two_mol_idx(int i, int j)
{
	/* absolute index in the two-molecule-pair storage 
	   i >= j                                         */

	if ( i < j )
	{
		puts("Warning! two_mol_idx Function Behavior Undefined...");
	}
	else
	{
		puts("2_mol_idx Good");
	}


	int abs_idx = 0; 
	for (int x = 0; x < i; x++)
	{
		abs_idx += x + 1;
	}
	abs_idx += j;
	return abs_idx;
}


int three_mol_idx(int i, int j, int k)
{
	/* absolute index in the three-molecule-pair storage 
	   i >= j >= k                                       */

	if ( !(i>=j && j>=k) )
	{
		puts("Warning! three_mol_idx Function Behavior Undefined...");
	}
	else
	{
		puts("3_mol_idx Good");
	}

	int abs_idx = 0; 
	for (int x = 0; x < i; x++)
	{
		for (int y = 0; y <= x; y++)
		{
			abs_idx += y + 1;
		}
	}

	for (int z = 0; z < j; z++)
	{
		abs_idx += z + 1;
	}
	abs_idx += k;
	return abs_idx;
}
'''