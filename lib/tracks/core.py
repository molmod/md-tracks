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

import numpy, os


__all__ = [
    "dump_track", "TrackNotFoundError", "load_track",
]


class Error(Exception):
    pass


class TrackNotFoundError(Error):
    pass


def complete_slice(sub):
    return slice(sub.start or 0, sub.stop or sys.maxint, sub.step or 1)

class TracksSplitter(object):
    """Efficiently convert any kind of trajectory file into separate tracks.

    After creating an object of this class, one calls the method dump_row
    for each record in the trajectory. When done, one calls the method finalize.
    """
    total_buffer_size = 10*1024*1024 # the size of the buffers, defaults to 1MB
    dot_interval = 50 # print a dot on screen, each 50 steps.

    def __init__(self, filenames):
        """Initialize the serialization procedure.

        Arguments:
            directory -- the directory where the tracks are stored
            filenames -- the file names of the individual tracks

        This function creates the directory and initializes the buffers for the
        serialization procedure.
        """
        for filename in filenames:
            directory = os.path.dirname(filename)
            if not os.path.exists(directory):
                os.mkdir(directory)
        self.filenames = filenames
        self.buffer_length = self.total_buffer_size/len(filenames)
        self.buffers = numpy.zeros((self.buffer_length, len(filenames)), float)
        self.counter = 0
        self.rowcount = 0
        self.flushcount = 0
        # cleanup existing stuff
        for filename in self.filenames:
            if os.path.isfile(filename):
                os.remove(filename)

    def dump_row(self, row):
        """Dump a record to the track files.

        The elements from row correspons to the file names in self.filenames.
        Each time the buffers are filled, they are written to disk and emptied.
        """
        self.rowcount += 1
        if self.rowcount % self.dot_interval == 0:
            log(".", False)
        self.buffers[self.counter] = row
        self.counter += 1
        if self.counter == self.buffer_length:
            self._flush_buffers()

    def finalize(self):
        """Finalize the tracks created by this TracksSplitter instance.

        First the buffers are written to disk and emptied. Consequently each
        track is transformed into a single continuous pickled numpy array.
        """
        self._flush_buffers()
        if self.flushcount > 1:
            self._serialize()
        else:
            log("No further serialization required.")
        log.finalize()

    def _flush_buffers(self):
        if self.counter == 0:
            return
        log(str(self.rowcount), False)
        for filename, b in zip(self.filenames, self.buffers.transpose()):
            f = file(filename, 'a')
            b[:self.counter].dump(f)
            f.close()
        self.counter = 0
        self.flushcount += 1

    def _serialize(self):
        for filename in self.filenames:
            f = file(filename, 'r')
            blocks = []
            try:
                while True:
                    blocks.append(numpy.load(f))
            except EOFError:
                pass
            f.close()
            f = file(filename, 'w')
            numpy.concatenate(blocks).dump(f)
            f.close()
            log("serialized %s" % filename)


class TracksJoiner(object):
    """Performs essentially the inverse operation of TracksSplitter, i.e. it
    converts a set of tracks into a single trajectory file.

    Instances of this class are initialized with a set of track objects, and
    behave like an iterator that yields a rows containing corresponding values
    from each track.
    """
    total_buffer_size = 10*1024*1024
    dot_interval = 50

    def __init__(self, filenames):
        """Initializes a TracksJoiner object."""
        import sha, time, os
        self.prefix = sha.new(str(os.getpid())+" "+str(time.time())).hexdigest()
        self.buffer_length = self.total_buffer_size/len(filenames)
        self.filenames = filenames

    def cleanup(self):
        import shutil
        if os.path.isdir(self.prefix):
            shutil.rmtree(self.prefix)

    def _in_slice(self, b_index, sub):
        row_start = b_index*self.buffer_length
        if row_start + self.buffer_length < sub.start: return False
        if row_start >= sub.stop: return False
        return True

    def yield_blocks(self, sub=slice(None)):
        """Yield blocks of rows of corresponding values from each track in self.filenames.

        First each track is de-serialized in small files stored into a temporary
        subdirectory. Then, one by one the temporary files are loaded into
        memory and the buffers from the TracksJoiner are reconstructed.
        """
        # complete the slice
        sub = complete_slice(sub)
        # a temporary working directory
        os.mkdir(self.prefix)
        # first de-serialize the track files
        try:
            self._de_serialize(sub)
            self.rowcount = 0
            for b_index in xrange(self.num_blocks):
                old_rowcount = self.rowcount
                # at which row does this block start:
                row_start = b_index*self.buffer_length
                # is this block included in the slice? If not, continue with the next loop
                if not self._in_slice(b_index, sub): continue
                # transform the total slice to a block slice
                block_start = sub.start - row_start
                if block_start < 0:
                    block_start += (-block_start/sub.step+1)*sub.step
                block_stop = sub.stop - row_start # always positive
                if self.num_blocks == 1:
                    filenames = self.filenames
                else:
                    filenames = ["%s/%i.%i" % (self.prefix, f_index, b_index) for f_index in xrange(len(self.filenames))]
                block = numpy.array([load_track(filename) for filename in self.filenames])
                block = block.transpose()
                block = block[block_start:block_stop:sub.step]
                yield block
                self.rowcount = old_rowcount + len(block)
                log(" %s " % self.rowcount, False)
        except:
            self.cleanup()
            log.finalize()
            raise
        self.cleanup()
        log.finalize()

    def yield_rows(self, sub=slice(None)):
        """Yield rows of corresponding values from each track in self.filenames."""
        for block in self.yield_blocks(sub):
            for row in block:
                yield row
                if self.rowcount % self.dot_interval == 0:
                    log(".", False)
                self.rowcount += 1

    def _de_serialize(self, sub=slice(None)):
        track_length = None
        for f_index, filename in enumerate(self.filenames):
            f = file(filename, 'r')
            t = numpy.load(f)
            f.close()
            if track_length is None:
                track_length = len(t)
            elif track_length != len(t):
                raise Error("Not all tracks are of the same length. (%s)" % filename)

            num_blocks = track_length/self.buffer_length+1
            if num_blocks == 1:
                log("De-serialization is not needed. The buffers are large enough.")
                break
            for b_index in xrange(num_blocks):
                if not self._in_slice(b_index, sub): continue
                block = t[b_index*self.buffer_length:(b_index+1)*self.buffer_length]
                block_filename = "%s/%i.%i" % (self.prefix, f_index, b_index)
                f = file(block_filename, 'w')
                block.dump(f)
                f.close()
            log("de-serialized %s" % filename)
        self.num_blocks = num_blocks


def dump_track(path, array):
    """Dump a numpy array into a track file, which is simply a pickled array."""
    array.dump(path)


def load_track(path):
    """Loads a track from file."""
    if not os.path.isfile(path):
        raise TrackNotFoundError("Track %s could not be found." % path)
    return numpy.load(path)


