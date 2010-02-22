# MD-Tracks is a statistical analysis toolkit for molecular dynamics
# and monte carlo simulations.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
# for Molecular Modeling (CMM), Ghent University, Ghent, Belgium; all rights
# reserved unless otherwise stated.
#
# This file is part of MD-Tracks.
#
# MD-Tracks is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# In addition to the regulations of the GNU General Public License,
# publications and communications based in parts on this program or on
# parts of this program are required to cite the following article:
#
# "MD-TRACKS: A productive solution for the advanced analysis of Molecular
# Dynamics and Monte Carlo simulations", Toon Verstraelen, Marc Van Houteghem,
# Veronique Van Speybroeck and Michel Waroquier, Journal of Chemical Information
# and Modeling, 48 (12), 2414-2424, 2008
# DOI:10.1021/ci800233y
#
# MD-Tracks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import sys


__all__ = [
    "Logger", "log", "usage_tail",
]


class Logger(object):
    """The Logger class regulates the output of the tracks scripts.

    It can be muted by setting log.verbose=False (see the log instance below).
    """
    def __init__(self, verbose=True, f=sys.stderr, old_newline=True):
        self.verbose = verbose
        self.f = f
        self.old_newline = old_newline

    def __call__(self, s, newline=True):
        if self.verbose:
            if newline:
                if not self.old_newline:
                    self.f.write("\n")
                self.f.write("%s\n" % s)
            else:
                self.f.write(s)
                self.f.flush()
            self.old_newline = newline

    def finish(self):
        if not self.old_newline and self.verbose:
            self.f.write("\n")

log = Logger()


usage_tail = """
The option -h prints out all available options

--------------------------------------------------------------------------------
TRACKS is developed at the Center for Molecular Modeling and is distributed
as free software, under the conditions of the GNU General public license,
version 3. For more information about TRACKS, visit:

http://molmod.ugent.be/code

If you like TRACKS, become a registered user via the website at no cost!
This will advance our funding proposals and facilitates the continuity of this
project.
--------------------------------------------------------------------------------
"""


