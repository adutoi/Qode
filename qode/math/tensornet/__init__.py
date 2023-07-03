#    (C) Copyright 2023 Anthony D. Dutoi
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

from .base     import evaluate, increment, extract, scalar_value
from .tensors  import primitive_tensor_wrapper, tensor_sum    # tensor_sum() can initialize an empty accumulator for += use
from .contract import contract                                # the only way to build a tensor_network
from .backends import dummy_backend, numpy_backend

dummy_tensor = primitive_tensor_wrapper(dummy_backend, copy_data=False)
np_tensor    = primitive_tensor_wrapper(numpy_backend, copy_data=False)
