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


from ccio.xyz import XYZReader, XYZWriter
from molmod.data import periodic
from molmod.units import angstrom, fs, parse_unit

import numpy, cPickle, os, glob, sys

__all__ = [
    "xyz_to_tracks", "cp2k_ener_to_tracks", "cpmd_traj_to_tracks",
    "tracks_to_xyz",
    "dump_track", "TrackNotFoundError", "load_track",
    "parse_slice", "parse_x_step", "parse_x_last", "parse_x_length"
    "dist_track", "bend_track", "dihed_track",
    "AtomFilter",
    "Logger", "log"
]


class Error(Exception):
    pass


class TracksSplitter(object):
    """Efficiently convert any kind of trajectory file into separate tracks.

    After creating an object of this class, one calls the method dump_row
    for each record in the trajectory. When done, one calls the method finalize.
    """
    total_buffer_size = 10*1024*1024 # the size of the buffers, defaults to 1MB
    dot_interval = 50 # print a dot on screen, each 50 steps.

    def __init__(self, directory, names):
        """Initialize the serialization procedure.

        Arguments:
            directory -- the directory where the tracks are stored
            filenames -- the file names of the individual tracks

        This function creates the directory and initializes the buffers for the
        serialization procedure.
        """
        if not os.path.exists(directory):
            os.mkdir(directory)
        self.filenames = [os.path.join(directory, name) for name in names]
        self.buffer_length = self.total_buffer_size/len(names)
        self.buffers = numpy.zeros((self.buffer_length, len(names)), float)
        self.counter = 0
        self.rowcount = 0
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
        self._serialize()

    def _flush_buffers(self):
        if self.counter == 0:
            return
        log(str(self.rowcount), False)
        for filename, b in zip(self.filenames, self.buffers.transpose()):
            f = file(filename, 'a')
            cPickle.dump(b[:self.counter], f)
            f.close()
        self.counter = 0

    def _serialize(self):
        for filename in self.filenames:
            f = file(filename, 'r')
            blocks = []
            try:
                while True:
                    blocks.append(cPickle.load(f))
            except EOFError:
                pass
            f.close()
            f = file(filename, 'w')
            cPickle.dump(numpy.concatenate(blocks), f)
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

    def next(self):
        return self

    def __iter__(self, sub=slice(None)):
        """Yield a row of corresponding values from each track in self.filenames.

        First each track is de-serialized in small files stored into a temporary
        subdirectory. Then, one by one the temporary files are loaded into
        memory and the buffers from the TracksJoiner are reconstructed. Given
        these buffers, the individual rows are reconstructed and yielded.
        """
        # interpret the slice
        start = sub.start or 0
        stop = sub.stop or sys.maxint
        step = sub.step or 1
        # a temporary working directory
        os.mkdir(self.prefix)
        try:
            # first de-serialize the track files
            self._de_serialize()
            counter = 0
            for b_index in xrange(self.num_blocks):
                blocks = [load_track("%s/%i.%i" % (self.prefix, f_index, b_index)) for f_index in xrange(len(self.filenames))]
                for row in zip(*blocks):
                    if counter >= start and counter < stop and (counter - start) % step == 0:
                        yield row
                        if counter % self.dot_interval == 0:
                            log(".", False)
                    counter += 1
                log(".", False)
        finally:
            import shutil
            shutil.rmtree(self.prefix)

    def _de_serialize(self):
        track_length = None
        for f_index, filename in enumerate(self.filenames):
            f = file(filename, 'r')
            t = cPickle.load(f)
            f.close()
            if track_length is None:
                track_length = len(t)
            elif track_length != len(t):
                raise Error("Not all tracks are of the same length. (%s)" % filename)

            num_blocks = len(t)/self.buffer_length+1
            for b_index in xrange(num_blocks):
                block = t[b_index*self.buffer_length:(b_index+1)*self.buffer_length]
                block_filename = "%s/%i.%i" % (self.prefix, f_index, b_index)
                f = file(block_filename, 'w')
                cPickle.dump(block, f)
                f.close()
            log("de-serialized %s" % filename)
        self.num_blocks = num_blocks


def xyz_to_tracks(filename, middle_word, destination, sub=slice(None), file_unit=angstrom, atom_indices=None):
    """Convert an xyz file into separate tracks."""
    xyz_reader = XYZReader(filename, sub, file_unit=file_unit)

    names = []
    if atom_indices is None:
        atom_indices = range(len(xyz_reader.numbers))
    else:
        atom_indices = list(atom_indices)
    for index in atom_indices:
        for cor in ["x", "y", "z"]:
            names.append("atom.%s.%07i.%s" % (middle_word, index, cor))

    tracks_splitter = TracksSplitter(destination, names)
    for title, coordinates in xyz_reader:
        tracks_splitter.dump_row(coordinates[atom_indices].ravel())
    tracks_splitter.finalize()


def cp2k_ener_to_tracks(filename, destination, sub=slice(None)):
    """Convert a cp2k energy file into separate tracks."""
    import itertools
    tracks_splitter = TracksSplitter(destination, ["step", "time", "kinetic_energy", "temperature", "potential_energy", "total_energy"])
    f = file(filename)
    for line in itertools.islice(f, sub.start, sub.stop, sub.step):
        row = [float(word) for word in line.split()[:6]]
        row[1] = row[1]*fs
        tracks_splitter.dump_row(row)
    f.close()
    tracks_splitter.finalize()


def cpmd_traj_to_tracks(filename, num_atoms, destination, sub=slice(None), atom_indices=None):
    """Convert a cpmd trajectory file into separate tracks.

    num_atoms must be the number of atoms in the system.
    """
    if atom_indices is None:
        atom_indices = range(num_atoms)
    else:
        atom_indices = list(atom_indices)
    tracks_splitter = TracksSplitter(destination, sum((
        ["atom.pos.%07i.x", "atom.pos.%07i.y", "atom.pos.%07i.z", "atom.vel.%07i.x", "atom.vel.%07i.y", "atom.vel.%07i.z"]
        for index in atom_indices
    ), []))

    f = file(filename)
    counter = 0
    row = []
    for line in f:
        words = line.split()[1:]
        for word in words:
            row.append(float(word))
        counter += 1
        if counter == num_atoms:
            tracks_splitter.dump_row(row)
            row = []
            counter = 0
    f.close()
    tracks_splitter.finalize()


def tracks_to_xyz(prefix, destination, symbols, sub=slice(None), file_unit=angstrom, atom_indices=None):
    """Converts a set of tracks into an xyz file."""
    if atom_indices is None:
        atom_indices = range(len(symbols))
    else:
        atom_indices = list(atom_indices)
    symbols = [symbols[index] for index in atom_indices]

    filenames = []
    for index in atom_indices:
        for c in 'xyz':
            filenames.append("%s.%07i.%s" % (prefix, index, c))

    f = file(destination, 'w')
    xyz_writer = XYZWriter(f, symbols, file_unit=file_unit)
    tracks_joiner = TracksJoiner(filenames)
    for row in tracks_joiner.yield_rows(sub):
        xyz_writer.dump("None", numpy.array(row).reshape((-1,3)))
    f.close()


def dump_track(path, array):
    """Dump a numpy array into a track file, which is simply a pickled array."""
    array.dump(path)


class TrackNotFoundError(Exception):
    pass


def load_track(path):
    """Loads a track from file."""
    if not os.path.isfile(path):
        raise TrackNotFoundError("Track %s could not be found." % path)
    return numpy.load(path)


def parse_slice(s):
    """Converts a text description of a slice into a slice object."""
    result = []
    for word, default in zip(s.split(":"), [0, sys.maxint, 1]):
        if word == '':
            result.append(default)
        else:
            result.append(float(word))
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

def parse_x_step(s, measure="time"):
    """Convert s into a discretization step.

    The argument s can be a track file that contains a equidistant x axis, or a
    distance between to subsequent data points multiplied by a unit, e.g. 1*fs.

    The return value is the discretization step in a.u.
    """
    def fn(x_axis):
        delta = x_axis[1:]-x_axis[:-1]
        if (delta[0] != delta).all():
            raise Error("The %s-axis is not equidistant. Is %s truly a %s-axis?" % (measure, s, measure))
        return delta[0]
    return _parse_x_track(s, fn)

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


def dist_track(prefix1, prefix2, sub):
    """Compute the distance between two atoms at each time step."""
    deltas = []
    for c in ['x','y','z']:
        first = load_track('%s.%s' % (prefix1, c))[sub]
        second = load_track('%s.%s' % (prefix2, c))[sub]
        deltas.append(second-first)
    distances = numpy.sqrt(sum(delta**2 for delta in deltas))
    return distances


def bend_track(prefix1, prefix2, prefix3, sub):
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
    dot[dot > 1] = 1
    dot[dot < -1] = -1
    angle = numpy.arccos(dot)
    return angle


def dihed_track(prefix1, prefix2, prefix3, prefix4, sub):
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
    dot[dot > 1] = 1
    dot[dot < -1] = -1
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

    def __init__(self, filter_atoms):
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

    def __call__(self, *test_indices):
        """Test wither one of the indexes belongs to the predefined set."""
        if self.filter_atoms is None:
            return True
        return len(self.filter_atoms.intersection(test_indices)) > 0


class Logger(object):
    """The Logger class regulates the output of the tracks scripts.

    It can be muted by setting log.verbose=False (see below).
    """
    def __init__(self, verbose=False, f=sys.stdout, old_newline=True):
        self.verbose = verbose
        self.f = f
        self.old_newline = old_newline

    def __call__(self, s, newline=True):
        if self.verbose:
            if newline:
                if not self.old_newline:
                    self.f.write("\n")
                self.f.write("%s\n" % s)
            else:
                self.f.write(s)
                self.f.flush()
            self.old_newline = newline

    def finalize(self):
        if self.old_newline and self.verbose:
            self.f.write("\n")

log = Logger()
