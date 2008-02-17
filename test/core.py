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


from common import *

from tracks.core import *
from tracks.log import log

import unittest, os, numpy


log.verbose = False


class TrackTestCase(BaseTestCase):
    def get_arrays(self):
        floats = [numpy.float32, numpy.float64]
        try:
            floats.append(numpy.float128)
        except AttributeError:
            pass
        try:
            floats.append(numpy.float96)
        except AttributeError:
            pass

        complexes = [numpy.complex64, numpy.complex128]
        try:
            complexes.append(numpy.complex256)
        except AttributeError:
            pass
        try:
            complexes.append(numpy.complex192)
        except AttributeError:
            pass

        ints = [numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64]

        # return a list of different types of arrays, one-dimensional and len=50
        return [
            numpy.random.normal(0,1,50).astype(f) for f in floats
        ] + [
            numpy.random.randint(0,10,50).astype(f) for f in ints
        ] + [
            (numpy.random.normal(0,1,50)+numpy.random.normal(0,1,50)*1.0j).astype(c) for c in complexes
        ]

    def test_load_dump(self):
        for rnd1 in self.get_arrays():
            dump_track("test", rnd1)
            rnd2 = load_track("test")
            self.assertArraysEqual(rnd1, rnd2)
            self.assertEqual(rnd1.dtype, rnd2.dtype)

    def test_append(self):
        for rnd1 in self.get_arrays():
            track = Track("test", clear=True)
            for index in xrange(10):
                track.append(rnd1[index*5:(index+1)*5])
            rnd2 = track.read()
            self.assertArraysEqual(rnd1, rnd2)
            self.assertEqual(rnd1.dtype, rnd2.dtype)

    def test_read_parts(self):
        for rnd1 in self.get_arrays():
            track = Track("test", clear=True)
            track.append(rnd1)
            rnd2 = []
            for index in xrange(10):
                rnd2.append(track.read(start=index*5, length=5))
            rnd2 = numpy.concatenate(rnd2)
            self.assertArraysEqual(rnd1, rnd2)
            self.assertEqual(rnd1.dtype, rnd2.dtype)

    def test_read_behind(self):
        for rnd1 in self.get_arrays():
            track = Track("test", clear=True)
            track.append(rnd1)
            rnd2 = track.read(50,10)
            self.assertEqual(len(rnd2), 0)
            self.assertEqual(rnd1.dtype, rnd2.dtype)


class MultiTrackTestCase(BaseTestCase):
    def get_buffers(self):
        return [
            numpy.random.normal(0,1,1000),
            numpy.random.normal(0,1,1000),
            numpy.random.randint(0,10,1000).astype(numpy.int32),
        ], ["test1", "test2", "test3"]

    def test_write(self):
        buffers, names = self.get_buffers()

        mtw = MultiTracksWriter(names, [b.dtype for b in buffers], buffer_size=5*1024)
        for row in zip(*buffers):
            mtw.dump_row(row)
        mtw.finalize()
        buffers_check = [load_track(name) for name in names]

        for b, b_check in zip(buffers, buffers_check):
            self.assertArraysEqual(b, b_check)
            self.assertEqual(b.dtype, b_check.dtype)

    def test_append(self):
        buffers, names = self.get_buffers()

        mtw = MultiTracksWriter(names, [b.dtype for b in buffers], buffer_size=5*1024)
        for row in zip(*buffers):
            mtw.dump_row(row)
        mtw.finalize()
        mtw = MultiTracksWriter(names, [b.dtype for b in buffers], buffer_size=5*1024, clear=False)
        for row in zip(*buffers):
            mtw.dump_row(row)
        mtw.finalize()

        buffers_check = [load_track(name) for name in names]

        for b, b_check in zip(buffers, buffers_check):
            b = numpy.concatenate([b,b])
            self.assertArraysEqual(b, b_check)
            self.assertEqual(b.dtype, b_check.dtype)

    def test_read(self):
        buffers, names = self.get_buffers()

        for b, name in zip(buffers, names):
            dump_track(name, b)
        mtr = MultiTracksReader(names, buffer_size=5*1024)
        buffers_check = list([] for index in xrange(len(buffers)))
        for row in mtr:
            for index, value in enumerate(row):
                buffers_check[index].append(value)
        buffers_check = [numpy.array(b, dtype) for b, dtype in zip(buffers_check, mtr.dtypes)]

        for b, b_check in zip(buffers, buffers_check):
            self.assertArraysEqual(b, b_check)
            self.assertEqual(b.dtype, b_check.dtype)




