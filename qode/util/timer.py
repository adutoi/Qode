#    (C) Copyright 2024 Anthony D. Dutoi
# 
#    This file is part of QodeApplications.
# 
#    QodeApplications is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    QodeApplications is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with QodeApplications.  If not, see <http://www.gnu.org/licenses/>.
#
import time


# A simple timer class that can accumulate time under several labels and
# report them as a percentage of the time elapsed since the timer was created.
class timer(object):
    def __init__(self):
        self._timings = {}
        self._t00 = time.time()        # the time of instance creation
    def start(self):
        self._t0 = time.time()         # start a timer
    def record(self, label):
        dt = time.time() - self._t0    # just the time since the timer last started
        if label not in self._timings:
            self._timings[label] = (0, 0.)
        n, t = self._timings[label]
        self._timings[label] = (n+1, t+dt)    # increment a category by the this interval
    def print(self, header=None):
        t_tot = time.time() - self._t00
        if header is not None:
            print(header)
        print("  Total time:  {:5.2f}".format(t_tot))
        t_accounted = 0
        for label,info in sorted(self._timings.items(), key=lambda item: item[1][1]):
            t_accounted += info[1]
            print("  {:30s}  {:5.2f}%  {:5d} calls".format(label, 100*info[1]/t_tot, info[0]))
        print("  {:5.2f}% accounted for".format(100*t_accounted/t_tot))
