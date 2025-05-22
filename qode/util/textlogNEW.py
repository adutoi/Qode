#    (C) Copyright 2024 Anthony D. Dutoi
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

# Some very simple (but powerful given the number of use cases below) functions and
# a class to facilitate indented printing and archiving of "screen" prints to a log.
# (The "screen" is the default behavior if the base function is the builtin print,
# but user can choose differently.)
#
# As backgound, the basic archiving pattern in our programming model is as follows:
#
# def some_workhorse_function(<application-specific arguments>, archive=None):
#     if archive is None:  archive = struct()
#     ...
#     archive.member = <some data>
#
# def some_workhorse_function(<application-specific arguments>, archive=None):
#     ...
#     if archive is not None:
#         archive.member = <some data>
#
# An archive is a struct that should be passed in from the outside, so that gathered
# data persists after an exception (if handled properly), for later debugging, mining,
# restarting, etc.  But even if a function is built with our archiving model in mind,
# they should be usable by codes that know nothing about it.  And such a mode can be
# also be used to suppress overly verbose archiving.  So there needs to be a default
# argument that means no archiving.  If the archive would be small, but with lots of
# calls to it, it can just make a garbage-collected dummy (first example).  Otherwise,
# be explicit throughout (second example).
#
# Meanwhile the generic template for a function that prints status updates is as follows:
#
# def some_workhorse_function(<application-specific arguments>, printout=print):
#     ...
#     printout("some timely message")
#
# This does not fit our archiving model, but it illustrates our general approach to
# printing, which has historically been conflated with archiving (sadly).  We do not
# demand or even request that functions print a text log.  All such printing is
# concieved of as timely status updates sent to the "screen" for the immediate benefit
# of the user.  Archiving data is another matter (see above), and the archive can also
# have a copy of this status update log (see below).  The last argument should have the
# same interface as print and is expected to format the output in the same way as well
# (more on this shortly).  What can change is destination (file, string, internet stream),
# or it might be suppressed via the no_print function, to keep down overly verbose logs
# (this does not count as reformatting because it is more like multiplication by zero).
# The exception to the rule about not changing formatting is that the lines to be printed
# might be indented by the printout handler, to keep track of function call depth,
# since it is anticipated to all come out in one continuous stream (this can also be
# hacked to divert only some updates).  However, if anything else about the formatting
# changes then indented text (which is processed by an indpendent utility) will not come
# out as expected, and there can be discrepencies (beyond indentation) with the logger
# copy (see below).  The indented function below returns a print-like function that gives
# indented output.
#   For the record, the logging module was looked into and has not been categorically
# rejected, but it seemed not really suited to following the progress of a scientific
# application.  Every(!) message would be proceeded by one of only a few prefixes (like
# INFO), which we would constantly have to make decisions about (which would likely never
# be consistent).  And there does not seem to be easy indentation facilities.  Also the
# workings and the documentation are huuuge!, whereas this is dead simple and good enough.
#
# def some_workhorse_function(<application-specific arguments>, archive=None, printout=print):
#     if archive is None:  archive = struct()
#     ...
#     archive.member = <some data>
#     printout("some timely message")
#
# def some_workhorse_function(<application-specific arguments>, archive=None, printout=print):
#     if archive is None:  archive = struct()
#     archive.textlog    = logger(printout)
#     ...
#     archive.member = <some data>
#     archive.textlog("some timely message")
#
# It is the perogative of the function whether it saves its printed status updates to its
# archive, but if it chooses to do so, an instance of the logger class will pass through what
# it is given to the specified printout function (to file/string/stream, perhaps indented,
# perhaps suppressed ... all per the above discussion) and save a copy in a local string.
# A logger instance can be convered to a string, yeilding this value.  Our model is that if
# such a log is saved to an archive it should have the name textlog.  If a function stores
# text logs, this will happen even if screen printing us supressed.  If the user disagrees
# with the function's decision to store text logs, they can be deleted from the outside upon
# completion.  Indentation of messages sent to a given logger instance should be handled 
# manually.  There will be no external indentation, regardless of how the printout argument
# handles its indentation; it would not make sense because this one is inherently nested
# within other archives.



# could be useful publicly, but also called below to facilitate indenting and keeping records

def str_print(*objects, sep=" ", end="\n", indent=""):
    """ has the same interface as builtin print but puts results into a string, with potential indentation"""
    string = sep.join([str(object) for object in objects]) + end
    lines = string.split("\n")
    lines = [indent+line if (len(line)>0) else line for line in lines]
    return "\n".join(lines)

# public functions that aid in printing to the "screen"

def no_print(*objects, sep=" ", end="\n", flush=False):
    """ can be used to suppress printing logs without suppressing archiving of messages """
    pass

def indented(printout, indent="  "):
    """ take an object or function with same signature as builtin print and wrap it so that prints result indented """
    def _indented(*objects, sep=" ", end="\n", flush=True):
        string = str_print(*objects, sep=sep, end=end, indent=indent)
        printout(string, end="", flush=flush)
    return _indented

# class that simultaneously print to the "screen" and keep a record of such

class logger(object):
    """ instance has call interface of builtin print, sending messages to printout (same interface), saving archival copy """
    def __init__(self, printout):
        self._messages = ""
        self._printout = printout
    def __call__(self, *objects, sep=" ", end="\n", flush=True):
        message = str_print(*objects, sep=sep, end=end)
        self._printout(message, end="", flush=flush)
        self._messages += message
    def __str__(self):
        return self._messages
