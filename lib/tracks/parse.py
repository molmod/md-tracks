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


from tracks.core import TrackNotFoundError, load_track
from molmod.units import parse_unit

import sys


__all__ = [
    "parse_slice", "get_delta", "parse_x_step", "parse_x_last",
    "parse_x_length",
]


def parse_slice(s):
    """Converts a text description of a slice into a slice object."""
    result = []
    for word, default in zip(s.split(":"), [0, sys.maxint, 1]):
        if word == '':
            result.append(default)
        else:
            result.append(int(word))
    return slice(*result)


def _parse_x_track(s, fn, convert=parse_unit):
    try:
        # first try to read the file
        x_axis = load_track(s)
        return fn(x_axis)
    except TrackNotFoundError:
        # then interpret s as a measure with units
        try:
            return convert(s)
        except ValueError:
            raise Error("Can not open file %s and can not interpret %s." % (s, s))


def get_delta(x_axis):
    delta = x_axis[1:] - x_axis[:-1]
    if (delta[0] != delta).all():
        raise Error("The %s-axis is not equidistant. Is %s truly a %s-axis?" % (measure, s, measure))
    return delta[0]


def parse_x_step(s, measure="time"):
    """Convert s into a discretization step.

    The argument s can be a track file that contains a equidistant x axis, or a
    distance between to subsequent data points multiplied by a unit, e.g. 1*fs.

    The return value is the discretization step in a.u.
    """
    return _parse_x_track(s, get_delta)


def parse_x_last(s):
    """Convert s into the total time or interval size.

    The argument can be a track file with the x-axis, or the interval size
    multiplied by a unit, e.g. 20*ps.

    The return value is the total size in a.u.
    """
    def fn(x_axis):
        return x_axis[-1] - x_axis[0]
    return _parse_x_track(s, fn)


def parse_x_length(s):
    """Convert s into an array length

    The argument can be a track file or an integer. The return value is the size
    of the array or the given integer.
    """
    def fn(x_axis):
        return len(x_axis)
    return _parse_x_track(s, fn, int)

