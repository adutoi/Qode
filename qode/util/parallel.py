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
import queue
import multiprocessing

# This source file is a "stub" which should eventually contain a host of methods for dynamic parallelization



class resources(object):
	""" just a hint of things to come.  this should be modified and passed down as resources are more restricted at lower layers """
	def __init__(self,n_cores):
		self.n_cores = n_cores



def _serial(function, inputs):
	""" circumvents some buggy multiprocessing behavior when it is not needed ... good for debugging too """
	return [function(*i) for i in inputs]

class _function_wrapper(object):
	def __init__(self,f):  self.f = f
	def __call__(self,x):  return self.f(*x)

def _simple_multiprocessing(function, inputs, n_workers):
	wrapped_function = _function_wrapper(function)
	worker_pool = multiprocessing.Pool(n_workers)
	outputs = worker_pool.map(wrapped_function,inputs)
	worker_pool.terminate()		# Needed because repeated calls to parallelize_task will open many processes and exhaust unix file descriptors in /proc
	return outputs

def parallelize_task(function,inputs,n_workers=None):
	if n_workers is None:  n_workers = multiprocessing.cpu_count() - 2
	if n_workers<=1:
		outputs = _serial(function,inputs)
	else:
		outputs =        _simple_multiprocessing(function,inputs,n_workers)
		#outputs = _experimental_multiprocessing(function,inputs,n_workers)
	return outputs



class _f(object):
	def __init__(self, result_q, aggregator, action, other_arguments):
		self._result_q        = result_q
		self._aggregator      = aggregator
		self._action          = action
		self._other_arguments = other_arguments
	def __call__(self, inp):
		self._aggregate_results( self._action(*( inp + self._other_arguments )) )
		return None
	def _aggregate_results(self, out):
		try:
			inc = self._result_q.get_nowait()
		except queue.Empty:
			self._result_q.put(out)
		else:
			out = self._aggregator(out,inc)
			self._aggregate_results(out)

class aggregate(object):
	def __init__(self, action, other_arguments, aggregator, num_cores):
		"""\
		action is the function you are interested in executing many times and aggregating the results of
			its call signature must be action(*(inp+other_arguments)),
			where inp is a tuple of those arguments that are variable (single arguments must be wrapped as 1-tuples)
			and other_arguments is a tuple of those arguments that stay the same, passed as the next constructor argument
		aggregator is a function of two outputs that aggregates them into one output of the same type, which is returned.
			If the output of action is not amenable to repeated unordered pairwise aggregation, then this is the wrong utility to use.
		num_cores is self explanatory

		Usage pattern should look something like (for func(x,y,z), where x will vary)
			A = parallel.aggregate(func, (y,z), value_adder, 10)
			inputs = [(x,) for x resulting from some other calculation]
			start_val = a (likely empty) object of the type returned by func 
			end_val = A.run(start_val,inputs)	inputs could be any iteratable/iterator object
		"""
		self._workers    = multiprocessing.Pool(num_cores)
		self._result_q   = multiprocessing.Manager().Queue()
		self._function   = _f(self._result_q, aggregator, action, other_arguments)
		self._aggregator = aggregator
	def run(self, aggregate, inputs):					# Have to give a starting point for aggregation otherwise _finalize can fail if no inputs (has happened to me)
		self._result_q.put(aggregate)
		self._workers.map(self._function, inputs, chunksize=1)		# Maybe there is a better idea than using map, to avoid making a discarded output object?
		aggregate = self._finalize()					# Not the same object as incoming aggregate!
		self._workers.terminate()					# Needed because repeated instantiations of aggregate will open many processes and exhaust unix file descriptors in /proc (must be called after finalize?)
		return aggregate
	def _finalize(self):
		aggregate = self._result_q.get_nowait()
		while True:
			try:                 increment = self._result_q.get_nowait()
			except queue.Empty:  break
			else:                aggregate = self._aggregator(aggregate,increment)
		return aggregate







# What is going on here is that multiprocessing must serialize (pickle) the inputs and the functions, and this has horrible problems with nested function
# definitions.  I tried to get around this in a couple ways, but neither worked in my case.  The former solution is copied from
# http://stackoverflow.com/questions/3288595/multiprocessing-how-to-use-pool-map-on-a-function-defined-in-a-class
# but apparently, even though it works for lambdas (which the above default does not), it did not work for my stuff.


"""\
def _do_work(f, in_q, out_q):
	while True:
		i,x = in_q.get()
		if i is None:  break
		out_q.put( (i,f(x)) )


def _experimental_multiprocessing(function, inputs, n_workers):

	inputs_q  = multiprocessing.Queue(1)
	outputs_q = multiprocessing.Queue()

	workers = [multiprocessing.Process(target=_do_work, args=(function, inputs_q, outputs_q)) for _ in range(n_workers)]

	for worker in workers:
		worker.daemon = True
		worker.start()

	for i,x in enumerate(inputs):  inputs_q.put( (i,x) )
	for   _ in range(n_workers):   inputs_q.put( (None,None) )
	results  = [ outputs_q.get() for _ in range(len(inputs)) ]

	for worker in workers:
		worker.join()

	outputs = [ f for i,f in sorted(results)]

	return outputs



if __name__ == '__main__':
	print(parallelize_task(lambda i: i * 2, [1, 2, 3, 4, 6, 7, 8],5))


"""


"""\

def _do_work(function_and_argument):
	f,x = function_and_argument
	return f(x)

def _experimental_multiprocessing(function, inputs, n_workers):
	worker_pool = multiprocessing.Pool(n_workers)
	outputs = worker_pool.map(_do_work,[(function,x) for x in inputs])
	worker_pool.terminate()
	return outputs
"""




