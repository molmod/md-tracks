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
from ccio.cp2k import CellReader
from molmod.units import angstrom, fs

import os, numpy, itertools


__all__ = [
    "xyz_to_tracks", "cp2k_ener_to_tracks", "cpmd_traj_to_tracks", "tracks_to_xyz",
]


def xyz_to_tracks(filename, middle_word, destination, sub=slice(None), file_unit=angstrom, atom_indexes=None, clear=True):
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

    mtw = MultiTracksWriter(filenames, clear=clear)
    for title, coordinates in xyz_reader:
        mtw.dump_row(coordinates[atom_indexes].ravel())
    mtw.finalize()


def cp2k_ener_to_tracks(filename, destination, sub=slice(None), clear=True):
    """Convert a cp2k energy file into separate tracks."""
    import itertools
    names = ["step", "time", "kinetic_energy", "temperature", "potential_energy", "total_energy"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtypes = [int, float, float, float, float, float]
    dtypes = [numpy.dtype(d) for d in dtypes]
    mtw = MultiTracksWriter(filenames, dtypes, clear=clear)
    f = file(filename)
    for line in itertools.islice(f, sub.start, sub.stop, sub.step):
        row = [float(word) for word in line.split()[:6]]
        row[1] = row[1]*fs
        mtw.dump_row(row)
    f.close()
    mtw.finalize()


def cpmd_ener_to_tracks(filename, destination, sub=slice(None), clear=True):
    """Convert a cp2k energy file into separate tracks."""
    import itertools
    names = ["step", "fict_kinectic_energy", "temperature", "potential_energy", "classical_energy", "hamiltonian_energy", "ms_displacement"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtypes = [int, float, float, float, float, float, float]
    dtypes = [numpy.dtype(d) for d in dtypes]
    mtw = MultiTracksWriter(filenames, dtypes, clear=clear)
    f = file(filename)
    for line in itertools.islice(f, sub.start, sub.stop, sub.step):
        row = [float(word) for word in line.split()[:7]]
        mtw.dump_row(row)
    f.close()
    mtw.finalize()

def cp2k_cell_to_tracks(filename, destination, sub=slice(None), clear=True):
    import itertools
    names = ["cell.a.x", "cell.a.y", "cell.a.z", "cell.b.x", "cell.b.y", "cell.b.z", "cell.c.x", "cell.c.y", "cell.c.z", "cell.a", "cell.b", "cell.c", "cell.alpha", "cell.beta", "cell.gamma"]
    filenames = list(os.path.join(destination, name) for name in names)
    mtw = MultiTracksWriter(filenames, clear=clear)
    cr = CellReader(filename)
    for cell in itertools.islice(cr, sub.start, sub.stop, sub.step):
        norms = numpy.sqrt((cell**2).sum(axis=0))
        alpha = numpy.arccos(numpy.clip(numpy.dot(cell[:,1],cell[:,2])/norms[1]/norms[2], -1,1))
        beta = numpy.arccos(numpy.clip(numpy.dot(cell[:,2],cell[:,0])/norms[2]/norms[0], -1,1))
        gamma = numpy.arccos(numpy.clip(numpy.dot(cell[:,0],cell[:,1])/norms[0]/norms[1], -1,1))
        mtw.dump_row(numpy.concatenate([cell.transpose().ravel(), norms, [alpha, beta, gamma]]))
    mtw.finalize()


def cpmd_traj_to_tracks(filename, num_atoms, destination, sub=slice(None), atom_indexes=None, clear=True):
    """Convert a cpmd trajectory file into separate tracks.

    num_atoms must be the number of atoms in the system.
    """
    if atom_indexes is None:
        atom_indexes = range(num_atoms)
    else:
        atom_indexes = list(atom_indexes)
    names = sum((
        ["atom.pos.%07i.x" % index, "atom.pos.%07i.y" % index, "atom.pos.%07i.z" % index, "atom.vel.%07i.x" % index, "atom.vel.%07i.y" % index, "atom.vel.%07i.z" % index]
        for index in atom_indexes
    ), [])
    filenames = list(os.path.join(destination, name) for name in names)
    mtw = MultiTracksWriter(filenames, clear=clear)

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
