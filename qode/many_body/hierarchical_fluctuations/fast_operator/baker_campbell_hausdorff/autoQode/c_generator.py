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


index_dict = { 'i':'No', 'j':'No', 'k':'No', 'l':'No', 'm':'No', 'n':'No', 'o':'No', \
	           'a':'Nv', 'b':'Nv', 'c':'Nv', 'd':'Nv', 'e':'Nv', 'f':'Nv', 'g':'Nv' }

def two_index_to_position( list_of_index ):
	if len(list_of_index) != 2:
		raise AssertionError
	#
	#
	piece2 = " + " + list_of_index[1]
	piece1 =         list_of_index[0] + " * " + index_dict[list_of_index[1]]
	return piece1 + piece2

def four_index_to_position( list_of_index ):
	if len(list_of_index) != 4:
		raise AssertionError

	#
	#  Assume this indexing scheme:
	#    Absolute Position in Array = idx0 * ( num_orb_1 * num_orb_2 * num_orb_3 ) + idx1 * ( num_orb_2 * num_orb_3 ) + idx2 * num_orb_3 + idx3 
	#  Generate this string...
	#
	piece4 = " + " + list_of_index[3]
	piece3 = " + " + list_of_index[2] + " * " + index_dict[list_of_index[3]]
	piece2 = " + " + list_of_index[1] + " * " + index_dict[list_of_index[2]] + " * " + index_dict[list_of_index[3]] 
	piece1 =         list_of_index[0] + " * " + index_dict[list_of_index[1]] + " * " + index_dict[list_of_index[2]] + " * " + index_dict[list_of_index[3]]
	# print(piece1 + piece2 + piece3 + piece4)
	return piece1 + piece2 + piece3 + piece4



def six_index_to_position( list_of_index ):
	if len(list_of_index) != 6:
		raise AssertionError

	piece6 = " + " + list_of_index[5]
	piece5 = " + " + list_of_index[4] + " * " + index_dict[list_of_index[5]]
	piece4 = " + " + list_of_index[3] + " * " + index_dict[list_of_index[4]] + " * " + index_dict[list_of_index[5]] 
	piece3 = " + " + list_of_index[2] + " * " + index_dict[list_of_index[3]] + " * " + index_dict[list_of_index[4]] + " * " + index_dict[list_of_index[5]]
	piece2 = " + " + list_of_index[1] + " * " + index_dict[list_of_index[2]] + " * " + index_dict[list_of_index[3]] + " * " + index_dict[list_of_index[4]] + " * " + index_dict[list_of_index[5]]
	piece1 =         list_of_index[0] + " * " + index_dict[list_of_index[1]] + " * " + index_dict[list_of_index[2]] + " * " + index_dict[list_of_index[3]] + " * " + index_dict[list_of_index[4]] + " * " + index_dict[list_of_index[5]]

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

def full_case_gen(index_sets, summand, func_num, prefix):
	from base import copy
	#
	for factor in summand.factors[1].factors:
		factor.bottom.paired_with = []		# Just making up
		factor.top.paired_with    = []		# a new member name.
	for factor in summand.factors[1].factors:
		factor.bottom.paired_with += [factor.top]
		factor.top.paired_with    += [factor.bottom]
	#
	# End hack
	op_indices = []
	for factor in summand.factors[1].factors:  
		op_indices += [factor.bottom, factor.top]
	char_bank = Sum._index_ordering(op_indices, index_sets, should_match=True)

	# print("CHAR BANK",char_bank)

	sum_indices = []
	loop_upper  = []
	loop_lower  = []
	for chars in char_bank:
		sum_indices += [  chars[0] ]
		loop_upper  += [ [chars[1]] ]
		loop_lower  += [ ['0'] ]
		for relate,other in chars[2:]:
			# if not flag:
			if relate == '<':
				loop_upper[-1] = [other]
			elif relate == '>':
				loop_lower[-1] = [other + '+1']
			# else:
			# 	if relate == '<':
			# 		loop_upper[-1] += [other]
			# 	elif relate == '>=':
			# 		loop_lower[-1] += [other]


	# print("TOP LEVEL SUM OBJ ========================")
	# print("INDICES =", sum_indices)
	# print("LOWERS  =", loop_lower)
	# print("UPPERS  =", loop_upper)

	outer_func_name = "void compute_contraction_" + str(func_num) + '('
	outer_func_body = '{\nomp_set_dynamic(0);\nomp_set_num_threads(nthd);\n#pragma omp parallel for\n'




	outer_func_body += "for ( int %s = %s; %s < %s; %s++ )\n{\n" %(sum_indices[0], loop_lower[0][0], sum_indices[0], loop_upper[0][0], sum_indices[0])


	for i in range(1, len(sum_indices)):
		# print(sum_indices[i],loop_lower[i], loop_upper[i])
		if len(loop_lower[i]) > 1:
			outer_func_body += 'int %s_begin = max( %d, ' %(sum_indices[i], len(loop_lower[i]))
			# temp_output = ''
			for item in loop_lower[i]:
				outer_func_body += " %s," %(item)
				# temp_output += " %s," %(item)
			outer_func_body = outer_func_body[:-1] + ' );\n'
			# outer_func_body += 'printf("' + '%d ' * len(loop_lower[i]) + '\\n", ' + temp_output[:-1] + ');\n'
		else:
			outer_func_body += 'int %s_begin = %s;\n' %(sum_indices[i], loop_lower[i][0])
		#
		if len(loop_upper[i]) > 1:
			outer_func_body += 'int %s_end   = min( %d, ' %(sum_indices[i], len(loop_upper[i])) 
			# temp_output = ''
			for item in loop_upper[i]:
				outer_func_body += " %s," %(item)
				# temp_output += " %s," %(item)
			outer_func_body = outer_func_body[:-1] + ' );\n'
			# outer_func_body += 'printf("' + '%d ' * len(loop_upper[i]) + '\\n", ' + temp_output[:-1] + ');\n'
		else:
			outer_func_body += 'int %s_end   = %s;\n' %(sum_indices[i], loop_upper[i][0])

		#
		if '[%s]' %(sum_indices[i]) in outer_func_body: ########################### DIRTY FIX
			return 0
		# outer_func_body += 'printf("'+ sum_indices[i] +'_begin = %d; ' + sum_indices[i] +'_end = %d\\n", ' + sum_indices[i] + '_begin, ' + sum_indices[i] + '_end);\n'
		outer_func_body += "for ( int %s = %s_begin; %s < %s_end; %s++ )\n{\n" %(sum_indices[i], sum_indices[i], sum_indices[i], sum_indices[i], sum_indices[i]) 



	target_block = summand.factors[1].code()
	# print("TARGET BLOCK =",target_block)

	target_block_code = ""
	target_indices    = []
	for item in target_block:
		target_block_code += item[0]
		target_indices    += item[1]
	# target_block_code += "_block"

	outer_func_body += "#pragma omp atomic\n"
	outer_func_body += target_block_code + "_block[" + index_to_position(target_indices) + "] += " 

	inner_sum = summand.factors[0]
	inner_index_sets, inner_summand = inner_sum._pass_through_helper(copy, [])

	X = inner_summand.factors[-2]
	T = inner_summand.factors[-1]

	for C in [X,T]:
		for pair in range(len(C.top)):
			C.bottom[pair].paired_with = []		# Just making up
			C.top[pair].paired_with    = []		# a new member name.
	for C in [X,T]:
		for pair in range(len(C.top)):
			C.bottom[pair].paired_with += [C.top[pair]]
			C.top[pair].paired_with    += [C.bottom[pair]]
	#
	# End hack
	coeff_indices = []					# all the index pairs of the Coeffs, with lowest order indices at the end
	for pair in range(max(len(X.top),len(T.top))):		# bottom and top of each have same length
		p = -1-pair
		if pair<len(T.top):     coeff_indices = [T.top[p]]    + coeff_indices
		if pair<len(X.top):     coeff_indices = [X.top[p]]    + coeff_indices
		if pair<len(T.bottom):  coeff_indices = [T.bottom[p]] + coeff_indices
		if pair<len(X.bottom):  coeff_indices = [X.bottom[p]] + coeff_indices
	inner_char_bank = Sum._index_ordering(coeff_indices, inner_index_sets)
	# print("INNER CHAR BANK =", inner_char_bank)
	inner_sum_indices = []
	inner_loop_upper  = []
	inner_loop_lower  = []


	for chars in inner_char_bank:
		inner_sum_indices += [ chars[0] ]
		inner_loop_upper  += [ [chars[1]] ]
		inner_loop_lower  += [ ['0'] ]
		for relate,other in chars[2:]:
			# if not flag:
			if relate == '<':
				inner_loop_upper[-1] = [other]
			elif relate == '>':
				inner_loop_lower[-1] = [other + '+1']
			# else:
			# 	if relate == '<':
			# 		inner_loop_upper[-1] += [other]
			# 	elif relate == '>=':
			# 		inner_loop_lower[-1] += [other]



	# print("INNER SUM OBJ ========================")
	# print("INDICES =", inner_sum_indices)
	# print("LOWERS  =", inner_loop_lower)
	# print("UPPERS  =", inner_loop_upper)
	#
	#
	inner_func_name = "double compute_inner_contraction_" + str(func_num) + '('
	inner_func_body = "{\ndouble contraction_value = 0.0;\n"

	# inner_func_body += "for ( int %s = %s; %s < %s; %s++ )\n{\n" %(inner_sum_indices[0], inner_loop_lower[0][0], inner_sum_indices[0], inner_loop_upper[0][0], inner_sum_indices[0])


	for i in range(len(inner_sum_indices)):
		# inner_func_body += "for ( int %s = %s; %s < %s; %s++ )\n" %(inner_sum_indices[i], inner_loop_lower[i], inner_sum_indices[i], inner_loop_upper[i], inner_sum_indices[i])
		if len(inner_loop_lower[i]) > 1:
			inner_func_body += 'int %s_begin = max( %d, ' %(inner_sum_indices[i], len(inner_loop_lower[i]))
			for item in inner_loop_lower[i]:
				inner_func_body += " %s," %(item)
			inner_func_body = inner_func_body[:-1] + ' );\n'
		else:
			inner_func_body += 'int %s_begin = %s;\n' %(inner_sum_indices[i],inner_loop_lower[i][0])
		#
		if len(inner_loop_upper[i]) > 1:
			inner_func_body += 'int %s_end   = min( %d, ' %(inner_sum_indices[i], len(inner_loop_upper[i])) 
			for item in inner_loop_upper[i]:
				inner_func_body += " %s," %(item)
			inner_func_body = inner_func_body[:-1] + ' );\n'
		else:
			inner_func_body += 'int %s_end   = %s;\n' %(inner_sum_indices[i], inner_loop_upper[i][0])
		#
		if '[%s]' %(inner_sum_indices[i]) in inner_func_body: ########################### DIRTY FIX
			return 0
		inner_func_body += "for ( int %s = %s_begin; %s < %s_end; %s++ )\n{\n" %(inner_sum_indices[i], inner_sum_indices[i], inner_sum_indices[i], inner_sum_indices[i], inner_sum_indices[i]) 


	inner_func_body += "contraction_value += "

	inner_terms = inner_summand
	for inner_term in inner_terms.factors:
		term_code = inner_term.code()
		if len(term_code[1]) > 0:
			inner_func_body += "%s[%s]*" %(term_code[0], index_to_position(term_code[1])) 
			inner_func_name += " double* %s," %(term_code[0])
			outer_func_name += " double* %s," %(term_code[0])
		else:
			inner_func_body += term_code[0] + '*'
	inner_func_body = inner_func_body[:-1] + ';\n'
	inner_func_body += '}\n' * (inner_func_body.count('{') -1 )
	inner_func_body += 'return contraction_value;\n}\n'

	for inner_index in sum_indices:
		inner_func_name += " int %s," %(inner_index)
	inner_func_name += ' int No, int Nv )'
		
	type_stripped = [ item for item in inner_func_name.split() if 'int' not in item and 'double' not in item ]
	# print(type_stripped)
	for item in type_stripped:
		if item.endswith('[][2]'):
			outer_func_body += item.split('[')[0] + ' '
		else:
			outer_func_body += item + ' '
	outer_func_body += ";\n"
	outer_func_body += '}\n' * outer_func_body.count('{') 

	outer_func_name +=  " int No, int Nv, double* " + target_block_code + "_block , int nthd )"

	# print(outer_func_body)
	# print(inner_func_body)

	# print("\n\n\n\n===============================================================================")
	# print("INNER FUNC NAME =", inner_func_name)
	# print("INNER FUNC BODY")
	# print(inner_func_body)

	# print("\n\n\n\n===============================================================================")
	# print("OUTER FUNC NAME =", outer_func_name)	
	# print("OUTER FUNC BODY")
	# print(outer_func_body)
	# print("===============================================================================")

	code_parts   = '\n' + inner_func_name + '\n'  + inner_func_body + '\n' + outer_func_name + '\n' + outer_func_body
	header_parts = inner_func_name + ';\n' + outer_func_name + ';\n'
	
	# print(code_parts)
	# print(header_parts)

	append_to_file(prefix+'.c', code_parts)
	append_to_file(prefix+'.h', header_parts)
	


def two_mult_gen(index_sets, summand, func_num, prefix):
	# raise AssertionError
	for factor in summand.factors[1].factors:
		factor.bottom.paired_with = []		# Just making up
		factor.top.paired_with    = []		# a new member name.
	for factor in summand.factors[1].factors:
		factor.bottom.paired_with += [factor.top]
		factor.top.paired_with    += [factor.bottom]
	#
	# End hack
	op_indices = []
	for factor in summand.factors[1].factors:  
		op_indices += [factor.bottom, factor.top]
	char_bank = Sum._index_ordering(op_indices, index_sets, should_match=True)


	sum_indices = []
	loop_upper  = []
	loop_lower  = []
	for chars in char_bank:
		sum_indices += [  chars[0] ]
		loop_upper  += [ [chars[1]] ]
		loop_lower  += [ ['0'] ]
		for relate,other in chars[2:]:
			# if not flag:
			if relate == '<':
				loop_upper[-1] = [other]
			elif relate == '>':
				loop_lower[-1] = [other + '+1']
			# else:
			# 	if relate == '<':
			# 		loop_upper[-1] += [other]
			# 	elif relate == '>=':
			# 		loop_lower[-1] += [other]

	# print("TOP LEVEL SUM OBJ ========================")
	# print("INDICES =", sum_indices)
	# print("LOWERS  =", loop_lower)
	# print("UPPERS  =", loop_upper)
	outer_func_name = "void compute_contraction_" + str(func_num) + '('
	outer_func_body = '{\nomp_set_dynamic(0);\nomp_set_num_threads(nthd);\n#pragma omp parallel for\n'

	
	# for i in range(len(sum_indices)):
	# 	outer_func_body += "for ( int %s = %s; %s < %s; %s++ )\n" %(sum_indices[i], loop_lower[i], sum_indices[i], loop_upper[i], sum_indices[i])
	
	outer_func_body += "for ( int %s = %s; %s < %s; %s++ )\n{\n" %(sum_indices[0], loop_lower[0][0], sum_indices[0], loop_upper[0][0], sum_indices[0])



	for i in range(1, len(sum_indices)):
		# print(sum_indices[i],loop_lower[i], loop_upper[i])
		if len(loop_lower[i]) > 1:
			outer_func_body += 'int %s_begin = max( %d, ' %(sum_indices[i], len(loop_lower[i]))
			# temp_output = ''
			for item in loop_lower[i]:
				outer_func_body += " %s," %(item)
				# temp_output += " %s," %(item)
			outer_func_body = outer_func_body[:-1] + ' );\n'
			# outer_func_body += 'printf("' + '%d ' * len(loop_lower[i]) + '\\n", ' + temp_output[:-1] + ');\n'
		else:
			outer_func_body += 'int %s_begin = %s;\n' %(sum_indices[i], loop_lower[i][0])

		#
		if len(loop_upper[i]) > 1:
			outer_func_body += 'int %s_end   = min( %d, ' %(sum_indices[i], len(loop_upper[i])) 
			# temp_output = ''
			for item in loop_upper[i]:
				outer_func_body += " %s," %(item)
				# temp_output += " %s," %(item)
			outer_func_body = outer_func_body[:-1] + ' );\n'
			# outer_func_body += 'printf("' + '%d ' * len(loop_upper[i]) + '\\n", ' + temp_output[:-1] + ');\n'
		else:
			outer_func_body += 'int %s_end   = %s;\n' %(sum_indices[i], loop_upper[i][0])
		#
		if '[%s]' %(sum_indices[i]) in outer_func_body: ########################### DIRTY FIX
			return 0
		outer_func_body += "for ( int %s = %s_begin; %s < %s_end; %s++ )\n{\n" %(sum_indices[i], sum_indices[i], sum_indices[i], sum_indices[i], sum_indices[i]) 


	mult_obj1, mult_obj2 = summand.factors
	target_block = mult_obj2.code()

	target_block_code = ""
	target_indices    = []
	for item in target_block:
		target_block_code += item[0]
		target_indices    += item[1]

	outer_func_body += "#pragma omp atomic\n"
	outer_func_body += target_block_code + "_block[" + index_to_position(target_indices) + "] += " 

	inner_terms = mult_obj1.factors
	for inner_term in inner_terms:
		term_code = inner_term.code()
		if len(term_code[1]) > 0:
			outer_func_body += "%s[%s]*" %(term_code[0], index_to_position(term_code[1])) 
			outer_func_name += " double* %s," %(term_code[0])
		else:
			outer_func_body += term_code[0] + '*'
	outer_func_body = outer_func_body[:-1] + ';\n'
	outer_func_body += '}\n' * outer_func_body.count('{')
	outer_func_name +=  " int No, int Nv, double* " + target_block_code + "_block , int nthd )"
	# print("TARGET BLOCK =", target_block_code, target_indices)
	# print("OUTER FUNC NAME =", outer_func_name)	
	# print("OUTER FUNC BODY")
	# print(outer_func_body)	
	code_parts   = '\n' + outer_func_name + '\n' + outer_func_body
	header_parts = outer_func_name + ';\n'

	append_to_file(prefix+'.c', code_parts)
	append_to_file(prefix+'.h', header_parts)
	
	# print(header_parts)
	# print(code_parts)

def three_coeff_gen(index_sets, summand, func_num, prefix):
	X = summand.factors[-2]
	T = summand.factors[-1]

	for C in [X,T]:
		for pair in range(len(C.top)):
			C.bottom[pair].paired_with = []		# Just making up
			C.top[pair].paired_with    = []		# a new member name.
	for C in [X,T]:
		for pair in range(len(C.top)):
			C.bottom[pair].paired_with += [C.top[pair]]
			C.top[pair].paired_with    += [C.bottom[pair]]

	coeff_indices = []					# all the index pairs of the Coeffs, with lowest order indices at the end
	for pair in range(max(len(X.top),len(T.top))):		# bottom and top of each have same length
		p = -1-pair
		if pair<len(T.top):     coeff_indices = [T.top[p]]    + coeff_indices
		if pair<len(X.top):     coeff_indices = [X.top[p]]    + coeff_indices
		if pair<len(T.bottom):  coeff_indices = [T.bottom[p]] + coeff_indices
		if pair<len(X.bottom):  coeff_indices = [X.bottom[p]] + coeff_indices
	char_bank = Sum._index_ordering(coeff_indices, index_sets)

	sum_indices = []
	loop_upper  = []
	loop_lower  = []
	for chars in char_bank:
		sum_indices += [  chars[0] ]
		loop_upper  += [ [chars[1]] ]
		loop_lower  += [ ['0'] ]
		for relate,other in chars[2:]:
			# if not flag:
			if relate == '<':
				loop_upper[-1] = [other]
			elif relate == '>':
				loop_lower[-1] = [other + '+1']
			# else:
			# 	if relate == '<':
			# 		loop_upper[-1] += [other]
			# 	elif relate == '>=':
			# 		loop_lower[-1] += [other]


	# print("TOP LEVEL SUM OBJ ========================")
	# print("INDICES =", sum_indices)
	# print("LOWERS  =", loop_lower)
	# print("UPPERS  =", loop_upper)
	outer_func_name = "void compute_contraction_" + str(func_num) + '('
	outer_func_body = '{\nomp_set_dynamic(0);\nomp_set_num_threads(nthd);\n#pragma omp parallel for\n'

	

	outer_func_body += "for ( int %s = %s; %s < %s; %s++ )\n{\n" %(sum_indices[0], loop_lower[0][0], sum_indices[0], loop_upper[0][0], sum_indices[0])

	for i in range(1, len(sum_indices)):
		# print(sum_indices[i],loop_lower[i], loop_upper[i])
		if len(loop_lower[i]) > 1:
			outer_func_body += 'int %s_begin = max( %d, ' %(sum_indices[i], len(loop_lower[i]))
			# temp_output = ''
			for item in loop_lower[i]:
				outer_func_body += " %s," %(item)
				# temp_output += " %s," %(item)
			outer_func_body = outer_func_body[:-1] + ' );\n'
			# outer_func_body += 'printf("' + '%d ' * len(loop_lower[i]) + '\\n", ' + temp_output[:-1] + ');\n'
		else:
			outer_func_body += 'int %s_begin = %s;\n' %(sum_indices[i], loop_lower[i][0])
		#
		if len(loop_upper[i]) > 1:
			outer_func_body += 'int %s_end   = min( %d,' %(sum_indices[i], len(loop_upper[i])) 
			# temp_output = ''
			for item in loop_upper[i]:
				outer_func_body += " %s," %(item)
				# temp_output += " %s," %(item)
			outer_func_body = outer_func_body[:-1] + ' );\n'
			# outer_func_body += 'printf("' + '%d ' * len(loop_upper[i]) + '\\n", ' + temp_output[:-1] + ');\n'
		else:
			outer_func_body += 'int %s_end   = %s;\n' %(sum_indices[i], loop_upper[i][0])
		#
		if '[%s]' %(sum_indices[i]) in outer_func_body:  ########################### DIRTY FIX
			return 0
		outer_func_body += "for ( int %s = %s_begin; %s < %s_end; %s++ )\n{\n" %(sum_indices[i], sum_indices[i], sum_indices[i], sum_indices[i], sum_indices[i]) 



	# for i in range(len(sum_indices)):
	# 	outer_func_body += "for ( int %s = %s; %s < %s; %s++ )\n" %(sum_indices[i], loop_lower[i], sum_indices[i], loop_upper[i], sum_indices[i])
	outer_func_body += "#pragma omp atomic\n"
	outer_func_body += 'K_block[0] += '

	inner_terms = summand.factors
	for inner_term in inner_terms:
		term_code = inner_term.code()
		if len(term_code[1]) > 0:
			outer_func_body += "%s[%s]*" %(term_code[0], index_to_position(term_code[1])) 
			outer_func_name += " double* %s," %(term_code[0])
		else:
			outer_func_body += term_code[0] + '*'
	outer_func_body = outer_func_body[:-1] + ';\n'
	outer_func_body += '}\n' * outer_func_body.count('{')
	outer_func_name +=  " int No, int Nv, double* K_block , int nthd )"

	# print("OUTER FUNC NAME =", outer_func_name)	
	# print("OUTER FUNC BODY")
	# print(outer_func_body)	
	code_parts   = '\n' + outer_func_name + '\n' + outer_func_body
	header_parts = outer_func_name + ';\n'

	append_to_file(prefix+'.c', code_parts)
	append_to_file(prefix+'.h', header_parts)
	
	# print(header_parts)
	# print(code_parts)


def build_shared_lib(file_main_name):
	#
	# USE THIS WHEN ALL CODE GENERATIONS ARE DONE!
	#
	import os
	os.system("gcc --std=c99 -fopenmp -c -fPIC " + file_main_name + ".c -I. -o " + file_main_name + ".o" )
	os.system("gcc -shared -fopenmp " + file_main_name  + ".o -o " + file_main_name + ".so")


def add_obj_call_code(add_obj, prefix):
	import os
	os.system("rm -f %s.c %s.h %s.o %s.so" %(prefix, prefix, prefix, prefix))
	append_to_file(prefix+'.c', '#include <stdio.h>\n#include "%s.h"\n#include "omp.h"\n' %(prefix))
	file_main_name = prefix
	print("gcc --std=c99 -fopenmp -c -fPIC " + file_main_name + ".c -I. -o " + file_main_name + ".o")
	print("gcc -shared -fopenmp " + file_main_name  + ".o -o " + file_main_name + ".so")
	ct = 0
	for item in add_obj.terms:
		# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		item.code(ct, prefix)
		ct += 1
	build_shared_lib(prefix)





