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


from tracks.core import MultiTracksReader, MultiTracksWriter
from ccio.xyz import XYZReader, XYZWriter
from molmod.units import angstrom, fs

import os, numpy, itertools


__all__ = [
    "xyz_to_tracks", "cp2k_ener_to_tracks", "cpmd_traj_to_tracks", "tracks_to_xyz",
]


def xyz_to_tracks(filename, middle_word, destination, sub=slice(None), file_unit=angstrom, atom_indexes=None):
    """Convert an xyz file into separate tracks."""
    xyz_reader = XYZReader(filename, sub, file_unit=file_unit)

    filenames = []
    if atom_indexes is None:
        atom_indexes = range(len(xyz_reader.numbers))
    else:
        atom_indexes = list(atom_indexes)
    for index in atom_indexes:
        for cor in ["x", "y", "z"]:
            filenames.append(os.path.join(destination, "atom.%s.%07i.%s" % (middle_word, index, cor)))

    mtw = MultiTracksWriter(filenames)
    for title, coordinates in xyz_reader:
        mtw.dump_row(coordinates[atom_indexes].ravel())
    mtw.finalize()


def cp2k_ener_to_tracks(filename, destination, sub=slice(None)):
    """Convert a cp2k energy file into separate tracks."""
    import itertools
    names = ["step", "time", "kinetic_energy", "temperature", "potential_energy", "total_energy"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtypes = [int, float, float, float, float, float]
    mtw = MultiTracksWriter(filenames)
    f = file(filename)
    for line in itertools.islice(f, sub.start, sub.stop, sub.step):
        row = [float(word) for word in line.split()[:6]]
        row[1] = row[1]*fs
        mtw.dump_row(row)
    f.close()
    mtw.finalize()


def cpmd_traj_to_tracks(filename, num_atoms, destination, sub=slice(None), atom_indexes=None):
    """Convert a cpmd trajectory file into separate tracks.

    num_atoms must be the number of atoms in the system.
    """
    if atom_indexes is None:
        atom_indexes = range(num_atoms)
    else:
        atom_indexes = list(atom_indexes)
    names = sum((
        ["atom.pos.%07i.x", "atom.pos.%07i.y", "atom.pos.%07i.z", "atom.vel.%07i.x", "atom.vel.%07i.y", "atom.vel.%07i.z"]
        for index in atom_indexes
    ), [])
    filenames = list(os.path.join(destination, name) for name in names)
    mtw = MultiTracksWriter(filenames)

    f = file(filename)
    counter = 0
    row = []
    for line in f:
        words = line.split()[1:]
        for word in words:
            row.append(float(word))
        counter += 1
        if counter == num_atoms:
            mtw.dump_row(row)
            row = []
            counter = 0
    f.close()
    mtw.finalize()


def tracks_to_xyz(prefix, destination, symbols, sub=slice(None), file_unit=angstrom, atom_indexes=None):
    """Converts a set of tracks into an xyz file."""
    if atom_indexes is None:
        atom_indexes = range(len(symbols))
    else:
        atom_indexes = list(atom_indexes)
    symbols = [symbols[index] for index in atom_indexes]

    filenames = []
    for index in atom_indexes:
        for c in 'xyz':
            filenames.append("%s.%07i.%s" % (prefix, index, c))

    f = file(destination, 'w')
    xyz_writer = XYZWriter(f, symbols, file_unit=file_unit)
    mtr = MultiTracksReader(filenames)
    for row in itertools.islice(mtr, sub.start, sub.stop, sub.step):
        xyz_writer.dump("None", numpy.array(row).reshape((-1,3)))
    f.close()
