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

from .base           import evaluate, increment, raw, scalar_value, shape, initialize_timer, print_timings, ContractionError
from .tensors        import backend_contract_path
from .tensors        import primitive_tensor as _primitive_tensor
from .tensors        import tensor_sum as _tensor_sum
from .backends       import dummy_backend, numpy_backend, tensorly_backend

print("""
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
U S I N G   E X P E R I M E N T A L   T E N S O R N E T   C O D E
""")

class primitive_tensor_factory(object):
    def __init__(self, backend, data_wrap=lambda d: d, data_unwrap=lambda d: d):
        self._backend     = backend
        self._data_wrap   = data_wrap
        self._data_unwrap = data_unwrap
    def init(self, data):
        raw_tensor = self._data_wrap(data)
        return _primitive_tensor(self._backend, raw_tensor)
    def zeros(self, shape=None):    # if shape is None, take it from first term added.
        return _tensor_sum(self._backend, shape=shape)
    def data(self, tensor):
        return self._data_unwrap(raw(tensor))

dummy_tensor = primitive_tensor_factory(dummy_backend.functions)
np_tensor    = primitive_tensor_factory(numpy_backend.functions)
tl_tensor    = primitive_tensor_factory(tensorly_backend.functions)
