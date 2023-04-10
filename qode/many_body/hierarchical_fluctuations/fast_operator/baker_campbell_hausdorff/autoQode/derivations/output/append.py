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

allpy = [ item[:-3] for item in os.listdir() if item.startswith('commute') and item.endswith('.py') ]


for item in allpy:
	print(item)
	words = item.split('_')
	# print(words[-2],words[-1])
	str_to_append = "\nfrom c_generator import add_obj_call_code\nadd_obj_call_code(expression, '%s')\n" %(words[-2] + '_' + words[-1])
	# print(str_to_append)
	f = open(item + '.py', 'r')
	original = f.read()
	f.close()

	original += str_to_append

	f = open(item+'.py','w')
	f.write(original)
	f.close()