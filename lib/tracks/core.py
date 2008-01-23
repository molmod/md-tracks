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


from tracks.log import log
from tracks import context

import numpy, os, struct


__all__ = [
    "Error", "TrackNotFoundError",
    "Track",
    "load_track", "dump_track",
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
        f = file(self.filename, "rb")
        header = f.read(self.header_size)
        f.close()
        if header[:8] != "TRACKS_1":
            raise Error("Wrong header: %s is not a correct track filename" % self.filename)
        return numpy.dtype(header[8:].replace("0",""))

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

    def read(self, start=0, length=-1):
        dtype, f = self._get_read_buffer(start)
        #return numpy.fromfile(f, dtype, length, '')
        # Don't use numpy.fromfile because it print annoying warning messages
        # on stderr when reading behind the end of the file.
        size = length*dtype.itemsize
        data = f.read(length*dtype.itemsize)
        if len(data) < size:
            length = len(data)/dtype.itemsize
        return numpy.ndarray(shape=(length,), dtype=dtype, buffer=data)

    def append(self, data):
        if len(data.shape) != 1:
            raise Error("Only 1-dimensional arrays can be stored in tracks.")
        f = self._get_append_buffer(data.dtype)
        data.tofile(f)


def load_track(filename):
    return Track(filename).read()


def dump_track(filename, data):
    Track(filename, clear=True).append(data)


class MultiTracksReader(object):
    def __init__(self, filenames, buffer_size=None, dot_interval=None):
        if buffer_size is None:
            buffer_size = context.default_buffer_size
        print buffer_size
        if dot_interval is None:
            dot_interval = context.default_dot_interval
        self.tracks = [Track(filename) for filename in filenames]
        self.dtypes = [track._get_header_dtype() for track in self.tracks]
        self.buffer_length = buffer_size/sum(dtype.itemsize for dtype in self.dtypes)
        self.dot_interval = dot_interval
        self.row_counter = 0

    def yield_buffers(self):
        buffer_counter = 0
        while True:
            start = buffer_counter*self.buffer_length
            log(" %i " % start, False)
            buffers = [track.read(start, self.buffer_length) for track in self.tracks]
            shortest = min(len(b) for b in buffers)
            longest = max(len(b) for b in buffers)
            if shortest == self.buffer_length:
                yield buffers
            else:
                buffers = [b[:shortest] for b in buffers]
                if longest != shortest:
                    log( " Not all tracks are of equal length! ", False)
                yield buffers
                break
            buffer_counter += 1
        end = buffer_counter*self.buffer_length + shortest
        log(" %i " % end, False)
        log.finalize()

    def yield_rows(self):
        for buffers in self.yield_buffers():
            for row in zip(*buffers):
                self.row_counter += 1
                if self.row_counter % self.dot_interval == 0:
                    log(".", False)
                yield row

    __iter__ = yield_rows


class MultiTracksWriter(object):
    def __init__(self, filenames, dtypes=None, buffer_size=None, dot_interval=None, clear=True):
        if buffer_size is None:
            buffer_size = context.default_buffer_size
        print buffer_size
        if dot_interval is None:
            dot_interval = context.default_dot_interval
        # make sure the files can be created
        for filename in filenames:
            directory = os.path.dirname(filename)
            if len(directory) > 0 and not os.path.exists(directory):
                os.makedirs(directory)
        if dtypes is None:
            self.dtypes = [numpy.dtype(float) for index in xrange(len(filenames))]
        else:
            if len(dtypes) != len(filenames):
                raise Error("len(dtypes) != len(filenames)")
            self.dtypes = dtypes
        self.tracks = [Track(filename, clear=clear) for filename in filenames]
        self.buffer_length = buffer_size/sum(dtype.itemsize for dtype in self.dtypes)
        self.buffers = [numpy.zeros(self.buffer_length, dtype) for dtype in self.dtypes]
        self.current_row = 0
        self.dot_interval = dot_interval
        self.row_counter = 0
        log(" 0 ", False)

    def _flush_buffers(self):
        if self.current_row == 0:
            return
        for b, t in zip(self.buffers, self.tracks):
            t.append(b[:self.current_row])
        log(" %i " % self.row_counter, False)
        self.current_row = 0

    def dump_row(self, row):
        if len(row) != len(self.buffers):
            raise Error("The row must contain len(self.buffers)=%i values." % len(self.buffers))
        for index, value in enumerate(row):
            self.buffers[index][self.current_row] = value
        self.current_row += 1
        self.row_counter += 1
        if self.row_counter % self.dot_interval == 0:
            log(".", False)
        if self.current_row == self.buffer_length:
            self._flush_buffers()

    def dump_buffers(self, buffers):
        self._flush_buffers()
        for b, t in zip(buffers, self.tracks):
            t.append(b)

    def finalize(self):
        self._flush_buffers()
        log.finalize()




