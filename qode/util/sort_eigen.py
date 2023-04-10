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
import numpy

def sort_eigen(np_out, order="ascending"):
    """ Takes the outputs of a numpy eigensolver and sorts them by the eigenvalue, without otherwise changing the format """
    if   order=="ascending"  or order=="<":  reverse = False
    elif order=="descending" or order==">":  reverse = True
    else:  raise ValueError("Unrecognized ordering requested")
    e, C = np_out
    e, C = zip(*sorted(zip( e, C.T.tolist() ), key=lambda pair: pair[0], reverse=reverse))
    return numpy.array(e), numpy.array(C).T


