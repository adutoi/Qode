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
from sys import stdout
from io import StringIO

"""\
A bit of philosophy.

Some of the capability of this system is not at all replaceable by simple output files.  In fact, this could even be used for checkpointing (I think).

The major difference between this and some global output file (in terms of what is printed), is that printing by this system is by default post-mortem.
This is ok, unless your calculation crashes and your output goes with it, or if you are anxious to know the status of your calculation.

My point of view on the latter is that the correct place for status updates is the screen, so that is why there is an easy way to have everything echoed
to the screen.  Anyone proficient in UNIX can easily re-pipe the output for easy checking remotely. However, that might still get tedious, and so I leave the
door open for more permanent configuration changes by allowing echoing of the log to a file.

The former complaint is more substantial, but I assert that the onus here is on the programmer to handle all exceptions in such a way that the output
structure is either printed or returned.  The alternative is to break locality and insist on file handles etc being handed down from above (which
could still be done with this system ... and nice and neatly)
"""



def no_print(*objects, sep=' ', end='\n'):  pass		# mimics interface tp builtin print()

def str_print(*objects, sep=' ', end='\n'):			# mimics interface tp builtin print()
	""" prints to a string and returns it """
	string_file = StringIO()
	print(*objects,sep=sep,end=end,file=string_file)
	string = string_file.getvalue()
	string_file.close()
	return string

def _indent_string(string,indentation):
	""" indents a string by a given amount """
	new_string = indentation
	if len(string)>0:
		for c in string[:-1]:
			if c=='\n':  new_string += '\n'+indentation
			else:        new_string +=   c
		new_string += string[-1]
	return new_string

def indent(*objects, sep=' ',indentation="  "):			# mimics interface tp builtin print()
	""" like str_print, but indents the result.  Since result will often be printed (adding \n by default), omit the end argument as this is easy to correct for on the outside """
	string = str_print(*objects,sep=sep,end='')
	string = _indent_string(string,indentation)
	return string

class indented(object):
	""" take an object or function with same signature as builtin print and wrap it so that prints result indented """
	def __init__(self,base_print,indentation="  "):
		self._base_print  = base_print
		self._indentation = indentation
	def __call__(self,*objects, sep=' ', end='\n'):		# mimics interface tp builtin print()
		self._base_print(indent(*objects,sep=sep),end=end)



class textlog(object):
	"""\
	 a call to this object has the same interface as builtin print, but it stores them as strings, keeping track of indentation and logs that are sublogs (but resolve to a string) 

	2. To keep a human-readable time-ordered log of what is happening. A textlog may be converted to a string
	   since its internal array contains only strings or other textlogs (so recursion bottoms out at strings only).  The nice thing about letting it contain
	   nested logs is that we have a built in way to adjust verbosity after the fact when the log is being explored.  In principle we can dump the rendered string
	   to a file when the calculation is finished, but in practice, we will probably want to pickle it for later exploration.
	   To facilitate real-time feed-back, echo file-streams may be given that dump the log in real time to the screen or another file.
	   The .sublog() method is intentionally simple.  The only difference is that the new textlog is linked into the
	   parent in its time-ordered location.  It is still a free-standing log of its own (and later may be identified to the parent as such by the unused "descriptor")

	The nice thing about this is how well it interfaces with code that knows nothing about it.  If the code can accept something with the signature of the builtin
	print, then a textlog can be passed in, or else maybe that code takes a file handle, in which case you pass in a StringIO and then use this to print it.
	On the other hand low-level functions may be written such that they can accept either a textlog or default to the builtin print, and the utilities above
	can still be used for indenting.  To get the most out of this though, some code will be written under the assumption that what comes in is a textlog, which
	can further be parsed out into sublogs.

	"""
	def __init__(self,echo=False,indentation="  ",indent_level=0):
		""" echo is a list of file-like streams where the log info is supposed to be printed to in real time """
		if   echo is False or echo is None:  echo = []
		elif echo is True:                   echo = [stdout]
		else:
			try:     echo = list(echo)
			except:  echo = [echo]
		self._echo = echo
		self._contents = []
		self._indentation = indentation
		self._indent_level = indent_level
	def __call__(self,*objects, sep=' ', end='\n', descriptor=None):		# mimics interface tp builtin print()
		indentation = self._indentation * self._indent_level
		string = str_print(*objects,sep=sep,end=end)
		string = _indent_string(string,indentation)
		self._contents += [(descriptor,string)]
		for out in self._echo:
			print(string,sep="",end="",file=out)
			out.flush()
	def __getitem__(self, index):
		return self._contents[index]
	def __iter__(self):
		return iter(self._contents)
	def __len__(self):
		return len(self._contents)
	def __reversed__(self):
		return reversed(self._contents)
	def __str__(self):
		string = ""
		for _,item in self._contents:  string += str(item)	# the discarded "descriptor" is in anticipation of generating outputs with dynamic levels of detail (that perhaps used above container methods)
		return string
	def _append(self,converts_to_string,descriptor=None):
		self._contents += [(descriptor, converts_to_string)]
	def indent(self):
		self._indent_level += 1
		return self
	def unindent(self):
		self._indent_level -= 1
		return self
	def sublog(self,indent=True):
		""" makes a new log but injected into present log """
		sub = textlog(self._echo,self._indentation,self._indent_level)
		self._append(sub)
		if indent:  sub.indent()
		return sub
	def write(self,filename):
		f = open(filename,"w")
		f.write(str(self)+'\n')
		f.close()

class output(object):
	"""\
	This class is intentionally very uncluttered. 
	1. To act as a container which is little more than a (nested) dictionary of data from a computation.  By being passed down from above and populated
	   as we go, it lowers the programmer overhead to keeping a rich record of everything that happened.  Once pickled, it can be conjured up and explored.
	   Therefore, it is absolutely intended that a user add or even delete members.  This is facilitated by the __call__ which can add a bunch of keyword members
	   in one line.

	This object is best suited to a return value.  A calling function needs to know nothing about what this object is, it behaves just like a dictionary from
	which to extract output.
	"""
	def __init__(self,**kwargs):
		self.__dict__.update(kwargs)
	def __call__(self,**kwargs):
		self.__dict__.update(kwargs)
		return self
