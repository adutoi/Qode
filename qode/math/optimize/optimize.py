#    (C) Copyright 2018, 2020, 2021 Anthony D. Dutoi and Yuhong Liu
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
import time
import numpy
from ...util.output import indent, indented, no_print, str_print



class no_DIIS_class(object):
	def __init__(self):                             pass
	def __call__(self, x, Dx, iteration, textlog):  return (x + Dx), (Dx@Dx)**(1/2)    # same signature as __call__ of DIIS_controller
	def print_info(self,textlog):                   pass
no_DIIS = no_DIIS_class()

class DIIS_controller(object):
	def __init__(self, start=0, max_dim=5, textlog=no_print):
		if max_dim<1:
			raise ValueError("DIIS Subspace Must Be Greater or Equal to One")
		if start<0:
			raise ValueError("DIIS Iteration Must At Least Start at First Cycle")
		# DIIS parameters
		self._start   = start		# iteration at which DIIS should take over (iterations start at 0)
		self._max_dim = max_dim		# size of the set of old vectors to store
		self._textlog = textlog		# If not no_print, sort of assuming it is identical to what was given to the optimizer
		# DIIS state
		self._D            = []		# List of prior step vectors at each prior position, as they would be without DIIS
		self._x            = []		# List of each prior position plus the unmodified step taken from that position
		self._subspace_dim = 0		# Dimensino of current subspace
		self._matrix       = []		# The error matrix (stored as lists because easy to work with, and matrix is small)
	def __call__(self, x, Dx, iteration, textlog):
		new_x,step_norm = no_DIIS(x, Dx, iteration, textlog)
		if iteration>=self._start:
			new_x = self._do_DIIS(new_x, Dx, textlog)	# Modify step using linear regression
			Dx = new_x - x
			step_norm = (Dx@Dx)**(1/2)
		return new_x,step_norm
	def print_info(self,textlog):
		textlog("DIIS will be used to accelerate convergence.")
		textlog(indent("DIIS will be invoked at Cycle", self._start+1))
		textlog(indent("DIIS Max subspace dimension =", self._max_dim))
	def release_vecs(self):						# Should be called by code that instantiates DIIS_controller, if desired (or just let fall out of scope)
		self._D      = None
		self._x      = None
		self._matrix = None
	def _do_DIIS(self, x, Dx, textlog):
		textlog("DIIS extrapolation ...")	# This already indented in log of caller ...
		if hasattr(self._textlog,"sublog"):  textlog = self._textlog.sublog()
		else:                                textlog = indented(self._textlog)
		textlog = indented(textlog)		# ... so indent this again

		self._subspace_dim += 1
		self._x += [  x ]
		self._D += [ Dx ]
		if self._subspace_dim>self._max_dim:
			self._subspace_dim -= 1
			self._shrink_arrays()				# Throws away earlier vectors and matrix columns and rows
		textlog("subspace dimension =", self._subspace_dim)

		self._increment_matrix()				# Inner product of new Dx with all the old ones
		textlog("residual metric:")
		textlog(indent(numpy.matrix(self._matrix)))
		weights = self._compute_weights()			# Upcast matrix to numpy, extend an invert to get weights
		textlog("linear regression weights:")
		textlog(indent(numpy.matrix([weights]).T))

		new_x = 0						# New position is linear combination of all old positions
		for w,x in zip(weights,self._x):  new_x = new_x + w*x	# A bit of extra overhead here until int+=space._member works, or I establish a null-vector object
		return new_x
	def _shrink_arrays(self):
		self._x.pop(0)
		self._D.pop(0)
		self._matrix.pop(0)
		for column in self._matrix:  column.pop(0)
	def _increment_matrix(self):
		last_row = [ ( self._D[-1] @ D ) for D in self._D ]
		self._matrix += [last_row]
		for row,last_element in zip(self._matrix[:-1],last_row[:-1]):  row += [last_element]
	def _compute_weights(self):
		extended_matrix        = -1.0 * numpy.matrix( numpy.ones(( self._subspace_dim+1, self._subspace_dim+1 )) )
		extended_matrix[-1,-1] =  0.0   # Last element set to zero
		for i in range(self._subspace_dim):
			for j in range(self._subspace_dim):
				extended_matrix[i,j] = self._matrix[i][j]
		inverse = numpy.linalg.inv(extended_matrix)
		weights = [ -inverse[i,-1] for i in range(self._subspace_dim) ]		# DO NOT INCLUDE THE LAST ELEMENT WHICH IS LAMBDA
		return weights



class Newton_step_class(object):
	def __init__(self):
		self.name = "Newton\'s method"
	def __call__(self, f, x0, textlog):
		textlog("Computing Newton step ...")
		textlog(indent("Computing gradient and inverse Hessian ..."))
		_, g0, invH0 = f(x0, value=False)		# The value (discarded), gradient and inverse of the Hessian
		textlog(indent("Applying inverse Hessian to gradient ..."))
		Dx = -invH0 @ g0				# The simple QN step, if DIIS not used
		return Dx
Newton_step = Newton_step_class()



def Optimize(f, x0, thresh, step, diis=no_DIIS, textlog=print):
	"""\
	x0 is the initial estimation of the location of the minimum or maximum of a scalar valued function represented by f; it is expected to be a member of a space.
	The algorithm converges when the step size, relative to the length of the vector x is less than thresh ... absolute scale is the responsibility of the caller.
	The function f has one mandatory argument, being x0 and displacements from it.  If all further optional arguments are not given it should return
	 a tuple which is, in order, the value of f at that point, the gradient, and the inverse of the Hessian.  The function should also take three optional
	 boolean arguments named, "value", "gradient", and "inv_Hessian", which, when False, means that the indicated return value is discarded and need not be
	 calculated.  Of course, there is nothing wrong with calculating it and returning it.  Furthermore, there is no reason not to return a cheap approximate
	 inverse Hessian, rather than the expensive real thing.  The gradient should be a member of a space, and the inverse of the Hessian should be an operator on
	 that space.
	"""
	textlog("==================================================================================")
	textlog("Optimization using {}.".format(step.name))
	textlog("Convergence when step-vector relative norm is less than {}.".format(thresh))
	diis.print_info(textlog)
	step_norm = float("inf")
	iteration = 0
	t_start = time.time()
	while step_norm>thresh:
		textlog("==================================================================================")
		textlog("Optimization Cycle:", iteration+1)			# For developers, iteration starts at 0, for users iteration starts at 1 (so print that)
		Dx = step(f, x0, indented(textlog))				# Compute the step from the base algorithm
		x1,step_norm = diis(x0,Dx,iteration,indented(textlog))		# Pass through the DIIS controller for potential modification
		t_cycle = time.time()
		x0 = x1								# Update ...
		iteration += 1							# ... for next iteration.
		step_norm /= (x0@x0)**(1/2)					# What fraction of present location was most recent step? (for convergence check)
		textlog("Relative step norm = {:.3e}".format(step_norm))	# Inform user how close we are getting
		textlog("Cycle Cumulative Time =", t_cycle - t_start)
	textlog("==================================================================================")
	textlog("Optimization Converged! ({} cycles)".format(iteration))
	textlog("==================================================================================")
	textlog("Computing function value at final point ...")
	v1, _, __ = f(x1, gradient=False, inv_Hessian=False)
	textlog("Cycle Cumulative Time =", time.time() - t_start)
	textlog("Final function value =", v1)
	textlog("==================================================================================")
	return (v1,x1)
