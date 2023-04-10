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
from qode.math import *


class _params1D(object):
     def __init__(self, origin, delta, min_index, n_points, dimension):
        self.origin    = origin
        self.delta     = delta
        self.min_index = min_index
        self.n_points  = n_points
        self.dimension = dimension


class _params2D(object):
    def __init__(self, x_origin, x_delta, x_min_index, xn_points, 
                       y_origin, y_delta, y_min_index, yn_points, dimension):
        self.x_origin      = x_origin
        self.x_delta       = x_delta
        self.y_origin      = y_origin
        self.y_delta       = y_delta
        self.x_min_index   = x_min_index
        self.y_min_index   = y_min_index
        self.xn_points     = xn_points
        self.yn_points     = yn_points
        self.dimension     = dimension

def _parse_domain(domain):
    if isinstance(domain,float):  return (0,domain), 2**8, 1
    elif isinstance(domain,tuple):
        xy_domain,n_points = domain
        if not isinstance(xy_domain[0], tuple):      # it expects a integer or float to work
            dimension = 1
            return xy_domain,n_points, dimension
        if isinstance(xy_domain[0], tuple):
            if len(xy_domain) == 2:
                dimension = 2
                return xy_domain, n_points, dimension
    else:
        domain,n_points = domain
        if isinstance(domain,float):  return (0,domain), n_points, 1
        else:                         return domain,n_points, 1   # the original input

def _direct_space_params(domain):
    domains,xy_n_points,dimension = _parse_domain(domain)
    if dimension == 1:
       n_points = xy_n_points
       origin,limit = domains
       delta = float(limit-origin) / n_points
       min_index = 0
       return _params1D(origin, delta, min_index, n_points, dimension)
    if dimension == 2:
       xn_points, yn_points = xy_n_points
       x_domain,y_domain  = domains 
       x_origin,x_limit = x_domain
       y_origin,y_limit = y_domain
       x_delta = float(x_limit-x_origin) / xn_points
       y_delta = float(y_limit-y_origin) / yn_points
       x_min_index = 0
       y_min_index = 0
       return _params2D(x_origin, x_delta, x_min_index, xn_points, 
                   y_origin, y_delta, y_min_index, yn_points, dimension)

def _reciprocal_space_params(domain,k_origin):
    kx_origin           = k_origin
    ky_origin           = k_origin
    domains,xy_n_points,dimension = _parse_domain(domain)
    if dimension == 1:
        n_points = xy_n_points
        origin,limit = domains
        delta = 1 / float(limit-origin)
        min_index = -n_points // 2
        return _params1D(k_origin, delta, min_index, n_points, dimension) 
    if dimension == 2:
        xn_points, yn_points = xy_n_points
        x_domain, y_domain   = domains
        x_origin,x_limit     = x_domain
        y_origin,y_limit     = y_domain
        x_delta = 1 / float(x_limit-x_origin)
        y_delta = 1 / float(y_limit-y_origin)
        x_min_index = -xn_points // 2		# integer division intended, works for even and odd  ...  make use of aliasing of python negative indices
        y_min_index = -yn_points // 2
        return _params2D(kx_origin, x_delta, x_min_index, xn_points, 
                   ky_origin, y_delta, y_min_index, yn_points, dimension)



def _populate_1D(func,origin,delta,min_index,n_points,data):
    # This is pretty good, but not 100% robust.  A dictionary is iterable but is not ordered.  Does not seem easy to separate those two without testing to see if func is a basestring instance?
    if   hasattr(func,'__call__') and not hasattr(func,'__iter__'):
        max_index = min_index + n_points
        print(func)
        for i in range(min_index,max_index):  data[i] = func( origin + i*delta )    # last point in domain aliased to first
    elif hasattr(func,'__iter__') and not hasattr(func,'__call__'):
        n = min(n_points,len(func))                     # is this *too* robust, or just convenient?
        for i in range(n):  data[i] = func[i]           # data should not include last point
    else:  raise error('No meaningful parse for function type given, or ambiguous (callable and iterable)')
    return data

def _populate_2D(func2D,x_origin,y_origin,x_delta,y_delta,x_min_index,y_min_index,xn_points,yn_points,data):
    if   hasattr(func2D,'__call__') and not hasattr(func2D,'__iter__'):
        x_max_index = x_min_index + xn_points
        y_max_index = y_min_index + yn_points
        #print(x_max_index, x_min_index, xn_points, x_origin, x_delta, "X")
        #print(y_max_index, y_min_index, yn_points, y_origin, y_delta, "Y")
        #print(len(data), "<== data")
        for i in range(x_min_index,x_max_index):
            for j in range(y_min_index,y_max_index):
                data[i,j] = func2D( x_origin + i*x_delta , y_origin + j*y_delta)
    elif hasattr(func2D,'__iter__') and not hasattr(func2D,'__call__'):
        # this will work if func2D is a list of lists, but what about 2D numpy arrays?  (if have to choose, go with list of lists ... this is supposed to be the general layer)
        nx = min(xn_points, len(func2D    ) )
        ny = min(yn_points, len(func2D[0] ) )
        # print(func2D)
        # print(type(func2D))
        for i in range(nx):
            for j in range(ny):
                data[i,j] = func2D[i][j]
    else:  raise error('No meaningful parse for function type given, or ambiguous (callable and iterable)')



# The field argument in these data types should not be taken seriously, just easier to type than numpy.float64 and numpy.complex128 and leaves door open for swapping out numpy if desired


class function(object):
    def __init__(self,domain,func=None,field=complex):
        self._domain = domain		# used as a compatibility marker later
        self.field = field
        p = _direct_space_params(domain)
        # Do this here to enforce data type consistency demanded by space code (will be cast to numpy's floating-point real/complex class internally)
        # ... and because func=None should give back an empty vector
        # In case of inconsistency, real data will be quietly upcast, which is what we usually would want, and a warning will be issued if the imaginary part is truncated.
        # The latter would be a pretty explicit blunder, giving a complex valued function or data to a this class when explicitly declared real.
        if p.dimension == 1:
            self._data = numpy.array([field(0)]*p.n_points)
            if func is not None:  _populate_1D(func,p.origin,p.delta,p.min_index,p.n_points,self._data)
        if p.dimension == 2:
            self._data = numpy.array( [[field(0)]*p.yn_points]*p.xn_points )
            if func is not None:  _populate_2D(func,p.x_origin,p.y_origin,p.x_delta,p.y_delta,p.x_min_index,p.y_min_index,p.xn_points,p.yn_points,self._data)



class _operator_base(object):          pass
class _recip_op_base(_operator_base):  pass
class _local_op_base(_operator_base):  pass

class direct_space_potential(_local_op_base):
    def __init__(self,domain,func,field=complex):
        self._domain = domain
        self.field = field
        p = _direct_space_params(domain)
        self.dimension = p.dimension
        if p.dimension == 1:
           self._data = numpy.array([field(0)]*p.n_points)
           _populate_1D(func,p.origin,p.delta,p.min_index,p.n_points,self._data)		# func=None is not allowed (no empty operators)
        if p.dimension == 2:
           #print(p.xn_points, "xn_points")
           #print(p.yn_points, "yn_points")
           self._data = numpy.array( [[field(0)]*p.yn_points]*p.xn_points )
           _populate_2D(func,p.x_origin,p.y_origin,p.x_delta,p.y_delta,p.x_min_index,p.y_min_index,p.xn_points,p.yn_points,self._data)

class direct_space_operator(_operator_base):
    def __init__(self,domain,func2D,field=complex):
        self._domain = domain
        self.field = field
        p = _direct_space_params(domain)
        self._data = numpy.array([[field(0)]*p.n_points]*p.n_points)
        _populate_2D(func2D,p.origin,p.delta,p.min_index,p.n_points,self._data)		# func=None is not allowed (no empty operators)

class reciprocal_space_potential(_recip_op_base,_local_op_base):	# dreaded diamond ok, b/c base classes all empty
    def __init__(self,domain,func,field=complex,k_origin=0):
        self._domain = domain
        self.field = field
        p = _reciprocal_space_params(domain,k_origin)
        if p.dimension == 1:
           self._data = numpy.array([field(0)]*p.n_points)
           _populate_1D(func,p.origin,p.delta,p.min_index,p.n_points,self._data)		# func=None is not allowed (no empty operators)
        if p.dimension == 2:
           self._data = numpy.array( [[field(0)]*p.yn_points]*p.xn_points )
           _populate_2D(func,p.x_origin,p.y_origin,p.x_delta,p.y_delta,p.x_min_index,p.y_min_index,p.xn_points,p.yn_points,self._data)

class reciprocal_space_operator(_recip_op_base):
    def __init__(self,domain,func2D,field=complex,k_origin=0):
        self._domain = domain
        self.field = field
        p = _reciprocal_space_params(domain,k_origin)
        self._data = numpy.array([[field(0)]*p.n_points]*p.n_points)
        _populate_2D(func2D,p.origin,p.delta,p.min_index,p.n_points,self._data)		# func=None is not allowed (no empty operators)



# For now, this is just a rip off of qode/math/numpy_space_traits.py, but the idea is to later generalize this so that the user can use other storage types.

"""\
The usual computer-science definition of the discrete FT is such that the high-index slots are best interpreted
physically as negative wavenumbers in a generalized two's-complement alias sense (where the last slot of a vector
is aliased to the index -1, etc.)  To enforce a symmetry between the two sides of the FT, we also demand that,
internally, the 0 index represents the origin of the position representation, and we treat the high-index
slots as negative positions, to the left of the origin (which makes no physical difference due to the periodicity).
The consquence of this is only a phase mask to be applied to the members of the wavenumber vector to represent
the true decomposition in the global coordinate frame.

The numpy DFTs are defined as:
Forward:  Gn =         sum[m=0 to N-1] exp(-i2pi * m*n/N) Fm
Inverse:  Fm = (1/N) * sum[n=0 to N-1] exp(+i2pi * m*n/N) Gn
Due to aliasing, the indices may be renamed such that 0 occurs in the middle of each summation.  If one
then also inserts the definition Xm = m*dX for some arbitrary dX, and then defines dK such that dX*dK = N
(physically motivated by the fact that if dX = L/N, then dK is 1/L), then the numpy
DFT may be rewritten as:
Forward:  Gn =         sum[Xm=-(N/2)dX to ((N-1)/2)*dX] exp(-i2pi * Xm*Kn) Fm
Inverse:  Fm = (1/N) * sum[Kn=-(N/2)dK to ((N-1)/2)*dK] exp(+i2pi * Xm*Kn) Gn
where integer division (round down) is assumed in the expressions to evaluate the limits (works for even and odd N).
Now making the association Fm = F(Xm) for some periodic function, the forward transform differs from the
discretized integral by a factor of dX, so we make the further association dX*Gn = G(Kn) to arrive at the physical FTs
Forward:  G(Kn) = sum[Xm=-(N/2)dX to ((N-1)/2)*dX] exp(-i2pi * Xm*Kn) F(Xm) * dX
Inverse:  F(Xm) = sum[Kn=-(N/2)dK to ((N-1)/2)*dK] exp(+i2pi * Xm*Kn) G(Kn) * dK
because we have recognized 1/(N*dX) = dK.
This then provides the recipe for converting between the numpy DFT and the physical FT.
If the input is a discretized function on a grid, the output of the numpy forward DFT needs to be scaled by dX to account
for the missing infinitesemal in the summation that replaces the integral.  If this output is then taken as input
for the inverse DFT, then one can either divide out the dX just multiplied and send this into the inverse DFT (which
then divides by N to recover the norm), or, equivalently but more conceptually, the result can be multiplied by dK
to account for the missing infinitesemal in that integral, and then the whole thing can be multiplied by N to 
counteract the division by N done in the numpy inverse DFT.

The last remaining tweaks are to consider that if the domain does not start at zero, then each element of the DFT is
off by a K-dependent phase (easy to correct) and that the entire wavefunction may be considered to have some implicit CoM motion
meaning that the K interval does not start at 0 ... for now, I'm skipping these concerns since I have no real use for the latter
and the former will cancel out in operators that perform the forward and then inverse transforms.
"""

class space_traits(object):
    def __init__(self,domain,field=complex):
        self._domain = domain		# used as a compatibility marker
        if field is float:  self.field = numpy.float64		# These are different than the fields stored in the data types due to numpy's casting
        else:               self.field = numpy.complex128
    def check_member(self,v):
        if not isinstance(v,function):  raise error("function instance expected")
        # all concerns about dimensionality should be taken care of by the next line?
        if v._domain!=self._domain:             raise error("cannot verify that domains are equivalent") 	# recommend reusing same domain object to get around float== test failing (or diff domains that parse same)
        p = _direct_space_params(self._domain)
        if p.dimension == 1:
            if type(v._data[0]) is not self.field:  raise error("vector contains scalars of wrong field type")	# Just a "spot check"? ... numpy insists array is homogeneous?
        if p.dimension == 2:
            if type(v._data[0][0]) is not self.field:  raise error("vector contains scalars of wrong field type")  # Just a "spot check"? ... numpy insists array is homogeneous?
    def check_lin_op(self,op):       
        if not isinstance(op,_operator_base):   raise error("_operator_base instance expected")
        # all concerns about dimensionality should be taken care of by the next line?
        if op._domain!=self._domain:            raise error("cannot verify that domains are equivalent") 	# recommend reusing same domain object to get around float== test failing
        if isinstance(op,_local_op_base):
            functionable = True
            p = _direct_space_params(self._domain)
            if p.dimension == 1:
                if type(op._data[0])   is not self.field:  raise error("diagonal array contains scalars of wrong field type")
            if p.dimension == 2:
                if type(op._data[0,0]) is not self.field:  raise error("operator matrix contains scalars of wrong field type")
        else:
            functionable = False
        return functionable
    @staticmethod
    def dot(v,w):
        p = _direct_space_params(v._domain)			# seems wierd to keep recomputing this.  store it in the vector?
        if p.dimension == 1:
            return p.delta * numpy.dot(v._data.conj(),w._data)
        if p.dimension == 2:
            total_integral = 0
            #print(len(v._data), len(w._data))
            for i in range(p.xn_points):
                for j in range(p.yn_points):
                    differential   =  p.x_delta * p.y_delta * v._data[i][j].conj() * w._data[i][j]
                    total_integral += differential
            return total_integral
    @staticmethod
    def add_to(v,w,n):
        if n==1:  v._data +=     w._data
        else:     v._data += n * w._data
    @staticmethod
    def scale(n,v):
        v._data *= n
    @staticmethod
    def copy(v):
        return function(domain=v._domain,func=v._data.copy(),field=v.field)
    @staticmethod
    def function_on_diags(func,op):
        func = numpy.vectorize(func)
        fop_data = func(op._data)
        if isinstance(op,_recip_op_base):
            return reciprocal_space_potential(domain=op._domain,func=fop_data,field=op.field)
        else:
            return direct_space_potential(domain=op._domain,func=fop_data,field=op.field)
    @staticmethod
    def act_on_vec(op,v):
        u_data = v._data
        p = _direct_space_params(v._domain)
        if p.dimension == 1:
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.fft(u_data)		# norm now off by 1/sqrt(N)
            if isinstance(op,_local_op_base):  u_data = numpy.multiply(op._data,u_data)
            else:                              u_data = numpy.dot(op._data,u_data)
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.ifft(u_data)		# ifft renormalizes automatically
            return function(domain=v._domain,func=u_data,field=v.field)
        if p.dimension == 2:
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.fft2(u_data)          # norm now off by 1/sqrt(N)
            if isinstance(op,_local_op_base):  u_data = numpy.multiply(op._data,u_data)
            else:                              u_data = numpy.dot(op._data,u_data)
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.ifft2(u_data)         # ifft renormalizes automatically
            return function(domain=v._domain,func=u_data,field=v.field)
    @staticmethod
    def back_act_on_vec(v,op):
        u_data = v._data
        p = _direct_space_params(v._domain)
        if p.dimension == 1:
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.fft(u_data)		# norm now off by 1/sqrt(N)
            if isinstance(op,_local_op_base):  u_data = numpy.multiply(op._data.conj(),u_data)
            else:                              u_data = numpy.dot(u_data.conj(),op._data).conj()
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.ifft(u_data)		# ifft renormalizes automatically
            return function(domain=v._domain,func=u_data,field=v.field)
        if p.dimension == 2:
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.fft2(u_data)          # norm now off by 1/sqrt(N)
            if isinstance(op,_local_op_base):  u_data = numpy.multiply(op._data.conj(),u_data)
            else:                              u_data = numpy.dot(u_data.conj(),op._data).conj()
            if isinstance(op,_recip_op_base):  u_data = numpy.fft.ifft2(u_data)         # ifft renormalizes automatically
            return function(domain=v._domain,func=u_data,field=v.field)
 



def id_operator(domain,field=complex):
    func = lambda x: 1
    return direct_space_potential(domain,func,field)

def position_operator(domain,field=complex):
    p = _direct_space_params(domain)
    if p.dimension == 1:
        func = lambda x: x
    if p.dimension == 2:
        func = lambda x,y: numpy.array([x,y])
    return direct_space_potential(domain,func,field)


# I am not generalizing this operator because the vector charater of the wavenumber operator in several dimensions, it conflicts with 
# the scalar character of the potential when _populate_2D fills out the grid/matrix. The code as it stands does not support operations 
# of vector operators like: p**2 = p * p 

def wavenumber_operator(domain,field=complex):                                                            
    """ In the most typical nomenclature this is the n operator, where n = k/2pi, with p = hbar*k  """    
    func = lambda n: n            									   
    return reciprocal_space_potential(domain,func,field)						   

def laplacian_operator(domain,field=complex):
    """ In the most typical nomenclature this is the n operator, where n = k/2pi, with p = hbar*k  """
    func = lambda nx,ny: (- (2*pi)**2 )*( (nx)**2 + (ny)**2 )
    return reciprocal_space_potential(domain,func,field)

def laplacian_operator_x2(domain,field=complex):
    """ In the most typical nomenclature this is the n operator, where n = k/2pi, with p = hbar*k  """
    func = lambda nx,ny: (- (2*pi)**2 )*(  (ny)**2 ) 
    return reciprocal_space_potential(domain,func,field)

def laplacian_operator_x1(domain,field=complex):
    """ In the most typical nomenclature this is the n operator, where n = k/2pi, with p = hbar*k  """
    func = lambda nx,ny: (- (2*pi)**2 )*(  (nx)**2 ) 
    return reciprocal_space_potential(domain,func,field)

def laplacian_operator_x1_x2(domain, m1, m2, field=complex):
    func = lambda nx,ny: (- (2*pi)**2 )*(  (nx)**2/m1 + (ny)**2/m2 )
    return reciprocal_space_potential(domain,func,field)

