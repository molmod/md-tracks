# Tracks provides tools for analyzing large trajectory files.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Tracks.
#
# Tracks is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Tracks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


from tracks.core import load_track

import numpy


__all__ = [
    "dist_track", "bend_track", "dihed_track",
    "AtomFilter",
]


def dist_track(prefix1, prefix2, sub=slice(None)):
    """Compute the distance between two atoms at each time step."""
    deltas = []
    for c in ['x','y','z']:
        first = load_track('%s.%s' % (prefix1, c))[sub]
        second = load_track('%s.%s' % (prefix2, c))[sub]
        deltas.append(second-first)
    distances = numpy.sqrt(sum(delta**2 for delta in deltas))
    return distances


def bend_track(prefix1, prefix2, prefix3, sub=slice(None)):
    """Compute the bending angle of three atoms at each time step."""
    deltas_a = []
    deltas_b = []
    for c in ['x','y','z']:
        first = load_track('%s.%s' % (prefix1, c))[sub]
        second = load_track('%s.%s' % (prefix2, c))[sub]
        third = load_track('%s.%s' % (prefix3, c))[sub]
        deltas_a.append(first-second)
        deltas_b.append(third-second)
    # compute the norms
    norm_a = numpy.sqrt(sum(delta_a**2 for delta_a in deltas_a))
    norm_b = numpy.sqrt(sum(delta_b**2 for delta_b in deltas_b))
    # normalize the vectors
    for delta_a in deltas_a:
        delta_a /= norm_a
    for delta_b in deltas_b:
        delta_b /= norm_b
    # calculate the dot product
    dot = sum(delta_a*delta_b for delta_a, delta_b in zip(deltas_a, deltas_b))
    dot = numpy.clip(dot, -1, 1)
    angle = numpy.arccos(dot)
    return angle


def dihed_track(prefix1, prefix2, prefix3, prefix4, sub=slice(None)):
    """Compute the dihedral angle of three atoms at each time step."""
    deltas_a = []
    deltas_b = []
    deltas_c = []
    for c in ['x','y','z']:
        first = load_track('%s.%s' % (prefix1, c))[sub]
        second = load_track('%s.%s' % (prefix2, c))[sub]
        third = load_track('%s.%s' % (prefix3, c))[sub]
        fourth = load_track('%s.%s' % (prefix4, c))[sub]
        deltas_a.append(first-second)
        deltas_b.append(third-second)
        deltas_c.append(fourth-third)
    # compute the norm of b
    norm_b = numpy.sqrt(sum(delta_b**2 for delta_b in deltas_b))
    # normalize the vectors b
    for delta_b in deltas_b:
        delta_b /= norm_b
    # project a and c on the plane orthogonal to b
    dot_ab = sum(delta_a*delta_b for delta_a, delta_b in zip(deltas_a, deltas_b))
    dot_cb = sum(delta_c*delta_b for delta_c, delta_b in zip(deltas_c, deltas_b))
    for delta_a, delta_b in zip(deltas_a, deltas_b):
        delta_a -= dot_ab*delta_b
    for delta_c, delta_b in zip(deltas_c, deltas_b):
        delta_c -= dot_cb*delta_b
    # compute the norms of a' and c'
    norm_a = numpy.sqrt(sum(delta_a**2 for delta_a in deltas_a))
    norm_c = numpy.sqrt(sum(delta_c**2 for delta_c in deltas_c))
    # normalize the vectors a' and c'
    for delta_a in deltas_a:
        delta_a /= norm_a
    for delta_c in deltas_c:
        delta_c /= norm_c
    # calculate the dot product and the angle
    dot = sum(delta_a*delta_c for delta_a, delta_c in zip(deltas_a, deltas_c))
    dot = numpy.clip(dot, -1, 1)
    angle = numpy.arccos(dot)
    swap = ((
        deltas_b[0] * (deltas_a[1] * deltas_c[2] - deltas_a[2] * deltas_c[1]) +
        deltas_b[1] * (deltas_a[2] * deltas_c[0] - deltas_a[0] * deltas_c[2]) +
        deltas_b[2] * (deltas_a[0] * deltas_c[1] - deltas_a[1] * deltas_c[0])
    ) > 0)*2-1
    angle *= swap
    return angle


class AtomFilter(object):
    """A tool to test whether some atoms belong to a user defined set."""

    def __init__(self, filter_atoms=None):
        """Initialize the atom filter.

        The argument filter_atoms can be a list of atom indexes or a string with
        comma-separated atom indexes.
        """
        if isinstance(filter_atoms, str):
            self.filter_atoms = frozenset(int(word) for word in filter_atoms.split(","))
        elif filter_atoms is None:
            self.filter_atoms = None
        else:
            self.filter_atoms = frozenset(filter_atoms)

    def __call__(self, *test_indexes):
        """Test wither one of the indexes belongs to the predefined set."""
        if self.filter_atoms is None:
            return True
        return len(self.filter_atoms.intersection(test_indexes)) > 0

