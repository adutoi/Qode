#    (C) Copyright 2018, 2023, 2025 Anthony D. Dutoi
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
import sys
import inspect
from . import struct

"""\
This module is dedicated to reading input out of python-syntax files that might contain only direct parameter
definitions or slightly more complex code aimed at defining parameters (there are no hardlimits on the code
complexity, but at a certain point it might make using this untenably brittle).

This has the dual benefits of making inputs quite flexible (especially for quick-and-dirty experiments concerning
how to specify parameters in an input deck) without having to write even a single line of parser code, but, at the
same time, the contents can be completely limited to parameter declarations (no need for import statements or other
"overhead").  All of this is accomplished while keeping the name of the input file command-line dynamic (which cannot
be accomplished easily using either imports of parameters from the main code, or imports of the main code from the
parameter code).

One of the more powerful aspects here is namespacing.  The author of the input file can directly import names
from any module it likes (including from the package that it is germane to), but also the author of the driver
code that accepts the input can have it excecuted in a custom namespace (which might just be its own namespace, 
passed via the globals() function).  This can be used to make available automatically names from other modules,
or perhaps specialized parser functions to lower the pain of writing input files.  Finally, a library of a common
set of such functions might eventually accompany this module, and those functions could be imported into the
executing namespace for convenience and uniformity.

The only drawback is that any parser errors will come back in "pythonese" which will not make sense to many
end users.  The proposed solution is another function (a plain text parser in the same module that defines input-level
namespaces?) that checks *simple* inputs for sanity (including whether they make physical sense, not just whether they
are parseable).  This check can be turned off by an intro line like (parse_check = False) that triggers the parse checker
to stop checking.  Then more advanced users can do things like embed loops, and other things that the parse checker 
does not like . . . try-except is also probably my friend in this battle.
"""



def from_string(string, write_to=None, namespace=None):
    """\
    Takes a string containing python code and executes it using the global namespace given as an argument
    (to expose names in modules like math, without them having to be referenced or imported in the string).
    The names defined inside the executed string are assigned to a namespace which is returned to the user.
    The identfiers in this space may be accessed using the . notation customary for builtin dictionaries.
    """
    if namespace is None:  namespace = {}    # if allow None through, the namespace of this module is used when exec is called
    if inspect.ismodule(namespace):
        namespace = namespace.__dict__       # in case it is a module of defintions for the parsing of inputs
    if write_to  is None:  write_to  = struct()
    exec(string, namespace, write_to)
    return write_to

def from_file(filename, write_to=None, namespace=None):
    """\
    Just a wrapper for read_input.from_string, but where the contents of an ASCII file are read as the string.
    """
    return from_string(open(filename).read(), write_to, namespace)

def from_argv(fields, write_to=None, namespace=None):
    """\
    Just a wrapper for read_input.from_string, but where a list of statements is executed (such as would be found in argv with no space around = signs)
    """
    string = ""
    for field in fields:  string += field + '\n'
    return from_string(string, write_to, namespace)

# This assumes that the command line is of the form
#    executable <params.in.py> [tokens interpreted as lines to append to <params.in.py>]
# If there are no defaults to define before parsing the parameters, then the simplest use is:
#     params = read_input.from_command_line()
# But if params is a pre-existing struct, either of the following will work
#     params.update(read_input.from_command_line())
#     read_input.from_command_line(write_to=params)
# All of these could be modified similarly to the first line becoming
#    params = read_input.from_command_line(namespace=parser_module)
# where parser_module contains defintions expected to be available in the input file environment
def from_command_line(write_to=None, namespace=None):
    params = from_file(sys.argv[1],  write_to=write_to, namespace=namespace)
    if len(sys.argv)>2:
        from_argv(sys.argv[2:], write_to=params, namespace=namespace)
    return params
