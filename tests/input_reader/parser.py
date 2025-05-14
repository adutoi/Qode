# This example module imports the some names and defines some
# functions that the user may want to have available automatically
# within the input file.  One might do this so that a module of
# parser functions can be imported to read the input without having
# to pollute the global namespace of the driver that calls the input
# reader.

from math import sqrt, pi, sin, cos

def string_to_matrix(string, convert=float, delimiter=None):
    # In the future, this should:
    #   1. by default ignore leading and trailing blank lines.  But the user should be
    #      able to turn this off, in case blank lines are meaningful.
    #   2. raise an exception if it encounters a reality counter to an indicated assumption,
    #      such as the "matrix" being square or even needing to be rectangular (the word
    #      "matrix" is anyway here already extended beyond needing to contain numbers).
    #      A handful of assumptions can be specified as booleans (already the convert argument
    #      can take care of a lot of this).
    lines = string.split("\n")
    matrix = []
    for line in lines:
        fields = line.split(delimiter)
        if convert is not None:
            fields = [convert(field) for field in fields]
        matrix += [fields]
    if len(matrix[ 0])==0:  matrix = matrix[1:]     # ignore first line if blank for input aestetics (user can always add another line to input)
    if len(matrix[-1])==0:  matrix = matrix[:-1]    # ignore last  line if blank for input aestetics (user can always add another line to input)
    return matrix
