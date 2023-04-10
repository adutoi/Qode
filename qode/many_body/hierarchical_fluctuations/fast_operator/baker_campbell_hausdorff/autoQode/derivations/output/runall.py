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


def wrapper(inputfile):
	os.system("python3.3 " + inputfile )

allpy = [ item for item in os.listdir() if item.startswith('commute') and item.endswith('.py') ]


# for item in allpy:
# 	wrapper(item)

myPool = Pool(20)
myPool.map(wrapper, allpy)
