# MD-Tracks is a statistical analysis toolkit for molecular dynamics
# and monte carlo simulations.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of MD-Tracks.
#
# MD-Tracks is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
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

from common import *

import tracks.vector as vector
import tracks.cell as cell

import numpy


__all__ = ["CellTestCase"]


class CellTestCase(BaseTestCase):
    def test_cell_basic(self):
        num = 10
        cell_track = cell.TrackCell(numpy.random.normal(0,1,(3,3,num)))
        cell_track.inv = numpy.array(cell_track.inv)
        for i in xrange(num):
            slab = cell_track.data[:,:,i]
            self.assertAlmostEqual(numpy.linalg.det(slab), cell_track.determinant[i])
            self.assertArraysAlmostEqual(numpy.linalg.inv(slab), cell_track.inv[:,:,i], 1e-5)

    def test_consistency1(self):
        num = 10
        cell_track = cell.TrackCell(numpy.random.normal(0,1,(3,3,num)))
        x = cell_track.to_fractional(vector.TrackVector(cell_track.data[:,0,:]))
        self.assertArrayAlmostConstant(x.data[0], 1, 1e-5)
        self.assertArrayAlmostZero(x.data[1], 1e-5)
        self.assertArrayAlmostZero(x.data[2], 1e-5)
        y = cell_track.to_fractional(vector.TrackVector(cell_track.data[:,1,:]))
        self.assertArrayAlmostZero(y.data[0], 1e-5)
        self.assertArrayAlmostConstant(y.data[1], 1, 1e-5)
        self.assertArrayAlmostZero(y.data[2], 1e-5)
        z = cell_track.to_fractional(vector.TrackVector(cell_track.data[:,2,:]))
        self.assertArrayAlmostZero(z.data[0], 1e-5)
        self.assertArrayAlmostZero(z.data[1], 1e-5)
        self.assertArrayAlmostConstant(z.data[2], 1, 1e-5)

    def test_consistency2(self):
        num = 10
        cell_track = cell.TrackCell(numpy.random.normal(0,1,(3,3,num)))
        r = vector.TrackVector(numpy.random.normal(0,1,(3,num)))
        k = cell_track.to_fractional(r)
        r_check = cell_track.from_fractional(k)
        self.assertArraysAlmostEqual(r.data, numpy.array(r_check.data), 1e-5)

    def test_shortest(self):
        num = 10
        cell_track = cell.TrackCell(numpy.random.normal(0,1,(3,3,num)))
        r = vector.TrackVector(numpy.random.normal(0,1,(3,num)))
        r_wrap = cell_track.shortest_vector(r)
        k_wrap = cell_track.to_fractional(r_wrap)
        self.assert_((numpy.array(k_wrap.data) <= 0.5).all())
        self.assert_((numpy.array(k_wrap.data) >= -0.5).all())



