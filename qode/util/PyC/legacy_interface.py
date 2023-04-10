import ctypes

def get_ctype(obj, data_type):
    if data_type == "int":
        return ctypes.c_longlong(obj)
    elif data_type == "int*":
        return obj.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong))
    elif data_type == "double*":
        return obj.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    else:
        raise AssertionError
