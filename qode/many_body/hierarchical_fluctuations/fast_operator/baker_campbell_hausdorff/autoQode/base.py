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
# advantage of external function is unbound name to give to _pass_through methods

def copy(indexed_base_object,index_replacements=[]):			# different signature than copy.copy !!!
	#print(indexed_base_object)
	return indexed_base_object._copy(index_replacements)

class indexed_base(object):
	""" basically, this base class is nothing more than documentation at this point """
	def __init__(self):                   pass
	def _copy(self, index_replacements):  raise NotImplementedError()
	def __str__(self):                    raise NotImplementedError()
