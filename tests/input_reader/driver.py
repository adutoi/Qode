#!/usr/bin/env python3
#
# Usage:
# ./driver.py params.in.py
# ./driver.py params.in.py 'number=(14,15)'   # the quotes would not be necessary if there were no parentheses

from sys import argv
from qode.util import read_input
import parser

# echo the input file for reference
print("============================")
print(open(argv[1]).read())
print("============================")

# extract the parameters in the input
# as an alternative to parser, could import define everything here and pass namespace=globals()
data = read_input.from_file(argv[1], namespace=parser)

# supplement or override from the command line
if len(argv)>2:
    read_input.from_argv(argv[2:], write_to=data)

# Here is what we got
print("number =", data.number)
print("letter =", data.letter)
print("irt2   =", data.irt2)
print("XY     =", data.XY)
print("norms  =", data.norms)
print()
