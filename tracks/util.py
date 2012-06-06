# -*- coding: utf-8 -*-
# MD-Tracks is a trajectory analysis toolkit for molecular dynamics
# and monte carlo simulations.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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
#--


import sys


__all__ = [
    "fix_slice", "AtomFilter",
]


def fix_slice(s):
    if s is None:
        return slice(0, sys.maxint, 1)
    else:
        return slice(s.start or 0, s.stop or sys.maxint, s.step or 1)


class AtomFilter(object):
    """A tool to test whether some atoms belong to a user defined set."""

    def __init__(self, filter_atoms=None):
        """Initialize the atom filter.

        The argument filter_atoms can be a list of atom indexes or a string with
        comma-separated atom indexes.
        """
        if isinstance(filter_atoms, str):
            if len(filter_atoms) == 0:
                self.filter_atoms = None
            else:
                self.filter_atoms = frozenset(int(word) for word in filter_atoms.split(","))
        elif filter_atoms is None:
            self.filter_atoms = None
        else:
            self.filter_atoms = frozenset(filter_atoms)

    def __call__(self, *test_indexes):
        """Test wether one of the indexes belongs to the predefined set."""
        if self.filter_atoms is None:
            return True
        return len(self.filter_atoms.intersection(test_indexes)) > 0


