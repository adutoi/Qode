# Comments are ignored, of course.

number = 13    # In one use case, this one will be overwritten on the command line.
letter = 't'

# Note that sqrt() is not defined by this file.  Since the driver code knows that
# this might be necessary, it passes in a namespace that defines this.
irt2  = 1/sqrt(2)



# Maybe you need to write some code to get the paremeters that your code wants in terms of some
# meta-parameters.  You can do anything you are used to doing in python (almost, see below).
theta = 15         # degrees ...
theta *= pi/180    # ... converted to radians.
XY = """
  3.0   9.0
 -1.0  -1.0
  2.5 -20.22
"""    # I like a string for this because it is easier than dealing with brackets and commas.
XY_ = []
for x,y in string_to_matrix(XY):       # The driver code has generously given a string_to_matrix function for us, ...
    x_ = x*cos(theta) + y*sin(theta)
    y_ = y*cos(theta) - x*sin(theta)
    XY_ += [(x_, y_)]
XY = XY_                               # ... and it only cares about the final contents of XY



# Maybe the calling code did not anticipate all of our needs; 
# we still can import what we need (perhaps even from the relevant package).
import numpy
norms = []
for vec in XY:
    norms += [numpy.linalg.norm(numpy.array(vec))]
# Not sure why this one-line version of the above fails.  Some limitation of exec()? ... still pretty powerful anyway
# norms = [numpy.linalg.norm(numpy.array(vec)) for vec in XY]
