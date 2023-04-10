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
from ...              import util
from ...math          import linear_inner_product_space
from ...math.optimize import Optimize, Newton_step, DIIS_controller, no_DIIS



class BCH_driver(object):
	def __init__(self, Tspace, BCH, H, invF, textlog):
		self.Tspace  = Tspace
		self.BCH     = BCH
		self.H       = H
		self.invF    = invF
		self.textlog = textlog
	def __call__(self, T_vec, value=True, gradient=True, inv_Hessian=True):
		textlog = self.textlog.sublog()						# Assume this is the textlog of the optimizer, so make a sublog for each call
		if inv_Hessian:  textlog.indent().indent()				# Match indentation inside optimizer loop
		#
		E_CCSD    = None							#
		Omega_vec = None							# Defaults if not requested
		invHess   = None							#
		#
		T = T_vec.v								# unwrap input
		#
		if value or gradient:
			textlog("Evaluating BCH expansion for  exp(-T) H exp(T)  ...")
			Omega = self.BCH.computeOmega(self.H, T, textlog.sublog())	# This is where all the time is spent!  (H, [H,T] , [[H,T],T] , ...)
			if gradient:  Omega_vec = self.Tspace.member(Omega)		# rewrap output
			#
			textlog("Projection:  E_CCSD = <0| [exp(-T) H exp(T)] |0>")
			E_CCSD  = self.BCH.Energy(Omega, textlog.sublog())		# So cheap at this point, do it regardless ...
			textlog('-------------------------')
			textlog("E_CCSD = {:.16e}".format(E_CCSD))
			textlog('-------------------------')
			if not value:  E_CCSD = None					# ... but censor it if not requested, for consistency/debugging
		#
		if inv_Hessian:
			textlog("Approximate inverse Hessian from perturbation theory  ~=  [1/(F-Eo)]")
			invHess = self.Tspace.lin_op(self.invF)
		#
		return E_CCSD, Omega_vec, invHess



def CCSD(H, T, invF, BCH, T_space_traits, textlog=util.textlog(echo=True), thresh=1e-15, diis_start=0, diis_subspace=5):
	if diis_start is None:  diis = no_DIIS
	else:                   diis = DIIS_controller(start=diis_start, max_dim=diis_subspace)
	Tspace = linear_inner_product_space(T_space_traits)
	E, Tvec = Optimize(f       = BCH_driver(Tspace, BCH, H, invF, textlog),
	                   x0      = Tspace.member(T),
	                   thresh  = thresh,
	                   step    = Newton_step,
	                   diis    = diis,
	                   textlog = textlog)
	return E, Tvec
