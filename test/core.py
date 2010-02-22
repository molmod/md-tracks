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


from common import *

from tracks.core import *
from tracks.log import log

import unittest, numpy


log.verbose = False


__all__ = ["TrackTestCase", "MultiTrackTestCase"]


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
                rnd2.append(track.read(slice(index*5, (index+1)*5)))
            rnd2 = numpy.concatenate(rnd2)
            self.assertArraysEqual(rnd1, rnd2)
            self.assertEqual(rnd1.dtype, rnd2.dtype)
            break

    def test_read_behind(self):
        for rnd1 in self.get_arrays():
            track = Track("test", clear=True)
            track.append(rnd1)
            rnd2 = track.read(slice(50,60))
            self.assertEqual(len(rnd2), 0)
            self.assertEqual(rnd1.dtype, rnd2.dtype)

    def test_read_into(self):
        sub = slice(10,30,2)
        for rnd1 in self.get_arrays():
            track = Track("test", clear=True)
            track.append(rnd1)
            destination = numpy.zeros(20,rnd1.dtype)
            track.read_into(destination, sub)
            self.assertArrayConstant(destination[10:],0)
            self.assertArraysEqual(destination[:10], rnd1[sub])


class MultiTrackTestCase(BaseTestCase):
    def get_data(self):
        dtype = numpy.dtype([("a", float, 2),("b", int, 1)])
        data = numpy.zeros(1000, dtype)
        data["a"] = numpy.random.normal(0,1,(1000,2))
        data["b"] = numpy.random.randint(0,10,1000)
        filenames = ["test1", "test2", "test3"]
        return data, filenames

    def dump_data(self, data, filenames):
        # this is the manual way of writing the data to disk
        counter = 0
        for name in data.dtype.names:
            sub_data = data[name]
            sub_dtype = data.dtype.fields[name][0]
            for flat_index in xrange(numpy.product(sub_dtype.shape,dtype=int)):
                index = numpy.unravel_index(flat_index, sub_dtype.shape)
                dump_track(filenames[counter], sub_data[(slice(None),)+index])
                counter += 1

    def read_data(self, dtype, size, filenames):
        # this is the manual way of reading the data from disk
        data = numpy.zeros(size, dtype)
        counter = 0
        for name in data.dtype.names:
            sub_data = data[name]
            sub_dtype = data.dtype.fields[name][0]
            for flat_index in xrange(numpy.product(sub_dtype.shape,dtype=int)):
                index = numpy.unravel_index(flat_index, sub_dtype.shape)
                tmp = load_track(filenames[counter])
                sub_data[(slice(None),)+index] = tmp
                counter += 1
        return data

    def compare_data(self, data1, data2):
        self.assertEqual(data1.dtype, data2.dtype)
        for name in data1.dtype.names:
            self.assertArraysEqual(data1[name], data2[name])

    def test_write(self):
        data, filenames = self.get_data()

        # write the data to disk, using the MultiTracksWriter
        mtw = MultiTracksWriter(filenames, data.dtype, buffer_size=5*1024)
        for row in data:
            mtw.dump_row(row)
        mtw.finish()

        # load it back manually
        data_check = self.read_data(data.dtype, len(data), filenames)

        # compare the original data with the data read from disk
        self.compare_data(data, data_check)

    def test_append(self):
        data, filenames = self.get_data()

        # write the data to disk, using the MultiTracksWriter
        mtw = MultiTracksWriter(filenames, data.dtype, buffer_size=5*1024)
        for row in data:
            mtw.dump_row(row)
        mtw.finish()
        mtw = MultiTracksWriter(filenames, data.dtype, buffer_size=5*1024, clear=False)
        for row in data:
            mtw.dump_row(row)
        mtw.finish()

        # load it back manually
        data_check = self.read_data(data.dtype, len(data)*2, filenames)

        # compare the original data with the data read from disk
        self.compare_data(data, data_check[:len(data)])
        self.compare_data(data, data_check[len(data):])

    def test_read(self):
        data, filenames = self.get_data()

        # write the data to disk, manually
        self.dump_data(data, filenames)

        # load the data back with the multi-tracks reader
        mtr = MultiTracksReader(filenames, data.dtype, buffer_size=5*1024)
        data_check = numpy.zeros(len(data), data.dtype)
        for index, row in enumerate(mtr):
            data_check[index] = row

        # compare the original data with the data read from disk
        self.compare_data(data, data_check)

    def test_read_sliced(self):
        data, filenames = self.get_data()
        sub = slice(10,120,13)

        # write the data to disk, manually
        self.dump_data(data, filenames)

        # load the data back with the multi-tracks reader
        mtr = MultiTracksReader(filenames, data.dtype, buffer_size=1024, sub=sub)
        data_check = numpy.zeros((sub.stop-sub.start-1)/sub.step+1, data.dtype)
        for index, row in enumerate(mtr):
            data_check[index] = row

        # compare the original data with the data read from disk
        self.compare_data(data[sub], data_check)


