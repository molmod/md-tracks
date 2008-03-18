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


from tracks.core import load_track
from tracks.util import fix_slice
from tracks.parse import parse_unit
from tracks.vector import TrackVector

from molmod.unit_cell import UnitCell

import numpy


__all__ = ["TrackCell"]


class TrackCell(object):
    @classmethod
    def from_cell_str(cls, unit_cell_str, sub=slice(None)):
        sub = fix_slice(sub)
        if "," in unit_cell_str:
            parameters = list(parse_unit(word) for word in unit_cell_str.split(",") if len(word) > 0)
            if len(parameters) == 1:
                a = parameters[0]
                single = numpy.array([[a,0,0],[0,a,0],[0,0,a]], float)
            elif len(parameters) == 3:
                a,b,c = parameters
                single = numpy.array([[[a,0,0],[0,b,0],[0,0,c]]], float)
            elif len(parameters) == 6:
                a,b,c,alpha,beta,gamma = parameters
                uc = UnitCell(
                    numpy.array([[1,0,0],[0,1,0],[0,0,1]], float),
                    numpy.array([True, True, True]),
                )
                uc.set_parameterst(numpy.array([a,b,c]), numpy.array([alpha,beta,gamma]))
                single = uc.cell
            elif len(parameters) == 9:
                single = numpy.array(parameters,float).reshape((3,3)).transpose()
            else:
                raise ValueError("If the --cell option contains comma's, one, three, six or nine value(s) are expected.")
            data = single.reshape((3,3,1))
        else:
            data = [
                [
                    load_track("%s.%s.%s" % (unit_cell_str, v, c), sub)
                    for v in "abc"
                ]
                for c in "xyz"
            ]
        return cls(data)

    def __init__(self, data):
        self.data = data
        # compute also the inverse:
        self.determinant = (
            data[0][0]*data[1][1]*data[2][2] + data[0][1]*data[1][2]*data[2][0] + data[1][0]*data[2][1]*data[0][2]
          - data[0][2]*data[1][1]*data[2][0] - data[1][0]*data[0][1]*data[2][2] - data[1][2]*data[2][1]*data[0][0]
        )
        self.inv = [
           [
              (data[1][1]*data[2][2]-data[1][2]*data[2][1])/self.determinant,
              (data[0][2]*data[2][1]-data[0][1]*data[2][2])/self.determinant,
              (data[0][1]*data[1][2]-data[0][2]*data[1][1])/self.determinant,
           ],
           [
              (data[1][2]*data[2][0]-data[1][0]*data[2][2])/self.determinant,
              (data[0][0]*data[2][2]-data[2][0]*data[0][2])/self.determinant,
              (data[0][2]*data[1][0]-data[0][0]*data[1][2])/self.determinant,
           ],
           [
              (data[1][0]*data[2][1]-data[1][1]*data[2][0])/self.determinant,
              (data[0][1]*data[2][0]-data[0][0]*data[2][1])/self.determinant,
              (data[0][0]*data[1][1]-data[0][1]*data[1][0])/self.determinant,
           ],
        ]

    def to_fractional(self, vector):
        return TrackVector([
            self.inv[0][0]*vector.data[0]+self.inv[0][1]*vector.data[1]+self.inv[0][2]*vector.data[2],
            self.inv[1][0]*vector.data[0]+self.inv[1][1]*vector.data[1]+self.inv[1][2]*vector.data[2],
            self.inv[2][0]*vector.data[0]+self.inv[2][1]*vector.data[1]+self.inv[2][2]*vector.data[2],
        ])

    def from_fractional(self, vector):
        return TrackVector([
            self.data[0][0]*vector.data[0]+self.data[0][1]*vector.data[1]+self.data[0][2]*vector.data[2],
            self.data[1][0]*vector.data[0]+self.data[1][1]*vector.data[1]+self.data[1][2]*vector.data[2],
            self.data[2][0]*vector.data[0]+self.data[2][1]*vector.data[1]+self.data[2][2]*vector.data[2],
        ])

    def shortest_vector(self, vector):
        frac = self.to_fractional(vector)
        for row in frac.data:
            row[:] -= row.round()
        return self.from_fractional(frac)

