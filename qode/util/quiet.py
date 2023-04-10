# This code adapted with gratitude from https://stackoverflow.com/questions/2828953/silence-the-stdout-of-a-function-in-python-without-trashing-sys-stdout-and-resto
# Since no copyright was claimed, we presume it to be in the public domain.
# The code was copied on 28.June.2018, and the user that posted it used the handle "Alex Martelli"

import sys
import contextlib

class DummyFile(object):
    def write(self, x):  pass

@contextlib.contextmanager
def quiet(send_to=DummyFile()):
    save = sys.stdout
    sys.stdout = send_to
    yield
    sys.stdout = save

# Usage:  given some function foo(), which prints to stdout without authorization, do
#
# with quiet():
#     foo()

# does not seem to work on psi4
# read comments below this answer on stackoverflow . . . seems this wont work if output coming from C.
# there is a way to do it and they are here:
# https://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python/17954769#17954769
