# MD-Tracks is a statistical analysis toolkit for molecular dynamics
# and monte carlo simulations.
# Copyright (C) 2007 - 2009 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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


from tracks.log import log
from tracks.util import fix_slice
from tracks import context

import numpy, os


__all__ = [
    "Error", "TrackNotFoundError",
    "Track",
    "load_track", "dump_track", "track_size",
    "MultiTracksReader", "MultiTracksWriter",
]


class Error(Exception):
    pass


class TrackNotFoundError(Error):
    pass


class Track(object):
    header_size = 14

    def __init__(self, filename, clear=False):
        self.filename = filename
        if clear:
            self.clear()

    def _init_buffer(self, dtype):
        f = file(self.filename, "wb")
        # write the header
        f.write("TRACKS_1") # file format and version
        f.write(dtype.str[:2]) # byte order and data type
        f.write("%04i" % dtype.itemsize) # the itemsize of the array in text format
        if f.tell() != self.header_size:
            raise Error("Inconsistent header size!")
        return f

    def _get_header_dtype(self):
        if not os.path.isfile(self.filename):
            raise TrackNotFoundError("File not found: %s" % self.filename)
        f = file(self.filename, "rb")
        header = f.read(self.header_size)
        f.close()
        if header[:8] != "TRACKS_1":
            raise Error("Wrong header: %s is not a correct track filename" % self.filename)
        return numpy.dtype(header[8:].replace("0",""))

    def _get_data_size(self):
        if not os.path.isfile(self.filename):
            raise TrackNotFoundError("File not found: %s" % self.filename)
        f = file(self.filename, "rb")
        f.seek(0, 2)
        size = f.tell() - self.header_size
        f.close()
        return size

    def _get_read_buffer(self, start):
        if not os.path.isfile(self.filename):
            raise TrackNotFoundError("File not found: %s" % self.filename)
        dtype = self._get_header_dtype()
        f = file(self.filename, "rb")
        f.seek(start*dtype.itemsize+14)
        return dtype, f

    def _get_append_buffer(self, dtype):
        if not os.path.isfile(self.filename):
            return self._init_buffer(dtype)
        else:
            dtype_file = self._get_header_dtype()
            if dtype != dtype_file:
                raise Error("The given data has dtype=%s, while the data in the track has dtype=%s" % (dtype, dtype_file))
            return file(self.filename, "ab")

    def clear(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def read(self, sub=None):
        sub = fix_slice(sub)
        dtype, f = self._get_read_buffer(sub.start)
        stop_bytes = min(sub.stop*dtype.itemsize, self._get_data_size())
        length = stop_bytes/dtype.itemsize - sub.start
        result = numpy.fromfile(f, dtype, length, '')[::sub.step]
        if not result.flags["C_CONTIGUOUS"]:
            result = result.copy()
        f.close()
        return result

    def read_into(self, destination, sub=None):
        sub = fix_slice(sub)
        dtype, f = self._get_read_buffer(sub.start)
        stop_bytes = min(sub.stop*dtype.itemsize, self._get_data_size())
        length = stop_bytes/dtype.itemsize - sub.start
        tmp = numpy.fromfile(f, dtype, length, '')[::sub.step]
        length = (length-1)/sub.step+1
        destination[:length] = tmp
        f.close()
        return length

    def append(self, data):
        if len(data.shape) != 1:
            raise Error("Only 1-dimensional arrays can be stored in tracks.")
        f = self._get_append_buffer(data.dtype)
        data.tofile(f)
        f.close()

    def size(self):
        dtype = self._get_header_dtype()
        size_bytes = self._get_data_size()
        return size_bytes/dtype.itemsize


def load_track(filename, sub=None):
    return Track(filename).read(sub)

def dump_track(filename, data):
    Track(filename, clear=True).append(data)

def track_size(filename):
    return Track(filename).size()


class MultiTrackBase(object):
    def init_buffer(self, buffer_size, dtype):
        # allocate the buffer array
        buffer_length = buffer_size/dtype.itemsize
        self.buffer = numpy.zeros(buffer_length, dtype)

    def init_tracks(self, filenames, dtype, clear=False):
        # create the tracks dictonary. it maps buffer array segments to filenames
        self.tracks = {}
        counter = 0
        for name in dtype.names:
            l = []
            self.tracks[name] = l
            sub_dtype = dtype.fields[name][0]
            for flat_index in xrange(numpy.product(sub_dtype.shape,dtype=int)):
                track = Track(filenames[counter], clear=clear)
                index = numpy.unravel_index(flat_index, sub_dtype.shape)
                l.append((index, track))
                counter += 1

    def _iter_fields(self, buffer=None):
        if buffer is None:
            buffer = self.buffer
        for name in self.buffer.dtype.names:
            sub_tracks = self.tracks[name]
            sub_buffer = buffer[name]
            for index, track in sub_tracks:
                column = sub_buffer[(slice(None),)+index]
                yield track, column


class MultiTracksReader(MultiTrackBase):
    def __init__(self, filenames, dtype, buffer_size=None, dot_interval=None, sub=slice(None)):
        MultiTrackBase.__init__(self)
        if buffer_size is None:
            buffer_size = context.default_buffer_size
        if dot_interval is None:
            dot_interval = context.default_dot_interval

        # some residual parameters
        self.dot_interval = dot_interval
        self.row_counter = 0
        self.sub = fix_slice(sub)

        self.init_buffer(buffer_size, dtype)
        self.init_tracks(filenames, dtype)
        self.init_shortest()

    def init_shortest(self):
        # compute the length of the shortest track in the reader
        self.shortest = None
        for tracks in self.tracks.itervalues():
            for i, track in tracks:
                size = track.size()
                if self.shortest is None or self.shortest > size:
                    self.shortest = size
        # take into account the slicing
        self.shortest = (min(self.shortest, self.sub.stop) - self.sub.start)/self.sub.step

    def iter_buffers(self):
        buffer_counter = 0
        while True:
            # determin the part that will be read from disk
            start = buffer_counter*len(self.buffer)*self.sub.step + self.sub.start
            if start >= self.sub.stop:
                break
            stop = min(start + len(self.buffer)*self.sub.step, self.sub.stop)

            # read the part slice(start, stop, step) from each track and store
            # it in the buffer array
            log(" %i " % start, False)
            first_size = None
            for track, column in self._iter_fields():
                size = track.read_into(column, slice(start, stop, self.sub.step))
                if first_size is None:
                    first_size = size
                elif first_size != size:
                    raise Error("Not all tracks are of equal length!")

            # yield the relevant part of the buffer array
            if size == len(self.buffer):
                yield self.buffer[:]
            else:
                yield self.buffer[:size]
                break
            buffer_counter += 1
        log(" %i " % stop, False)
        log.finish()

    def iter_rows(self):
        for buffer in self.iter_buffers():
            for row in buffer:
                self.row_counter += 1
                if self.row_counter % self.dot_interval == 0:
                    log(".", False)
                yield row

    __iter__ = iter_rows


class MultiTracksWriter(MultiTrackBase):
    def __init__(self, filenames, dtype, buffer_size=None, dot_interval=None, clear=True):
        MultiTrackBase.__init__(self)
        if buffer_size is None:
            buffer_size = context.default_buffer_size
        if dot_interval is None:
            dot_interval = context.default_dot_interval

        # make sure the files can be created
        for filename in filenames:
            directory = os.path.dirname(filename)
            if len(directory) > 0 and not os.path.exists(directory):
                os.makedirs(directory)

        self.init_buffer(buffer_size, dtype)
        self.init_tracks(filenames, dtype, clear)

        # some residual parameters
        self.current_row = 0
        self.dot_interval = dot_interval
        self.row_counter = 0
        log(" 0 ", False)

    def _flush_buffer(self):
        if self.current_row == 0:
            return
        for track, column in self._iter_fields():
            track.append(column[:self.current_row])
        log(" %i " % self.row_counter, False)
        self.current_row = 0

    def dump_row(self, row):
        #if row.dtype != self.buffer.dtype:
        #    raise Error("The row must have the same dtype as the internal buffer.")
        self.buffer[self.current_row] = row
        self.current_row += 1
        self.row_counter += 1
        if self.row_counter % self.dot_interval == 0:
            log(".", False)
        if self.current_row == len(self.buffer):
            self._flush_buffer()

    def dump_buffer(self, buffer):
        #if buffer.dtype != self.buffer.dtype:
        #    raise Error("The given buffer must have the same dtype as the internal buffer.")
        self._flush_buffer()
        for track, column in self._iter_fields(buffer):
            track.append(column)

    def finish(self):
        self._flush_buffer()
        log.finish()


