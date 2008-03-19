# Tracks provides tools for analyzing large trajectory files.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from tracks.core import TrackNotFoundError, Track, MultiTracksReader
from tracks.util import fix_slice

from molmod.units import parse_unit
from molmod.unit_cell import UnitCell

import sys, numpy


__all__ = [
    "Error", "parse_slice", "get_delta", "parse_x_step",
    "parse_x_duration", "parse_x_length", "yield_unit_cells",
]


class Error(Exception):
    pass


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
        x_track = Track(s)
        return fn(x_track)
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
    def fn(x_track):
        return get_delta(x_track.read(slice(10)))
    return _parse_x_track(s, fn)


def parse_x_duration(s):
    """Convert s into the total duration or interval size.

    The argument can be a track file with the x-axis, or the interval size
    multiplied by a unit, e.g. 20*ps.

    The return value is the total size in a.u.
    """
    def fn(x_track):
        size = x_track.size()
        return x_track.read(slice(size-1,size))[0] - x_track.read(slice(1))[0]
    return _parse_x_track(s, fn)


def parse_x_length(s):
    """Convert s into an array length

    The argument can be a track file or an integer. The return value is the size
    of the array or the given integer.
    """
    def fn(x_track):
        return x_track.size()
    return _parse_x_track(s, fn, int)


def yield_unit_cells(unit_cell_str, sub=None):
    sub = fix_slice(sub)
    if len(unit_cell_str) == 0:
        uc = UnitCell(
            numpy.array([[1,0,0],[0,1,0],[0,0,1]], float),
            numpy.array([False,False,False]),
        )
        while True:
            yield uc
    if "," in unit_cell_str:
        parameters = list(parse_unit(word) for word in unit_cell_str.split(",") if len(word) > 0)
        if len(parameters) == 1:
            a= parameters[0]
            uc = UnitCell(
                numpy.array([[a,0,0],[0,a,0],[0,0,a]], float),
                numpy.array([True, True, True]),
            )
        elif len(parameters) == 3:
            a,b,c = parameters
            uc = UnitCell(
                numpy.array([[a,0,0],[0,b,0],[0,0,c]], float),
                numpy.array([True, True, True]),
            )
        elif len(parameters) == 6:
            a,b,c,alpha,beta,gamma = parameters
            uc = UnitCell(
                numpy.array([[1,0,0],[0,1,0],[0,0,1]], float),
                numpy.array([True, True, True]),
            )
            uc.set_parameterst(numpy.array([a,b,c]), numpy.array([alpha,beta,gamma]))
        elif len(parameters) == 9:
            uc = UnitCell(
                numpy.array(parameters, float).reshape((3,3)),
                numpy.array([True, True, True]),
            )
        else:
            raise ValueError("If the --cell option contains comma's, one, three, six or nine value(s) are expected.")
        while True:
            yield uc
    else:
        filenames = ["%s.%s" % (unit_cell_str, suffix) for suffix in ["a.x", "a.y", "a.z", "b.x", "b.y", "b.z", "c.x", "c.y", "c.z"]]
        mtr = MultiTracksReader(filenames, sub=sub)
        uc = UnitCell(
            numpy.array([[1,0,0],[0,1,0],[0,0,1]], float),
            numpy.array([True, True, True]),
        )
        flat = uc.cell.ravel()
        for row in mtr:
            flat[:] = row
            uc.update_reciproke()
            yield uc



