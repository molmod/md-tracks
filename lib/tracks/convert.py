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


from tracks.core import MultiTracksReader, MultiTracksWriter
from molmod.io.xyz import XYZReader, XYZWriter
from molmod.io.cp2k import CellReader
from molmod.io.atrj import ATRJReader
from molmod.units import angstrom, fs

import os, numpy, itertools


__all__ = [
    "xyz_to_tracks", "cp2k_ener_to_tracks", "cpmd_traj_to_tracks", "tracks_to_xyz",
]


def yield_real_lines(f):
    for line in f:
        line = line[:line.find("#")]
        if len(line.strip()) > 0:
            yield line


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

    shape = (len(atom_indexes),3)
    dtype = numpy.dtype([("cor", float, shape)])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    for title, coordinates in xyz_reader:
        mtw.dump_row((coordinates[atom_indexes],))
    mtw.finalize()


def cp2k_ener_to_tracks(filename, destination, sub=slice(None), clear=True):
    """Convert a cp2k energy file into separate tracks."""
    names = ["step", "time", "kinetic_energy", "temperature", "potential_energy", "conserved_quantity"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtypes = [int, float, float, float, float, float]
    dtype = numpy.dtype([  (name, t, 1) for name, t in zip(names, dtypes)  ])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    f = file(filename)
    for line in itertools.islice(yield_real_lines(f), sub.start, sub.stop, sub.step):
        row = [float(word) for word in line.split()[:6]]
        row[1] = row[1]*fs
        mtw.dump_row(tuple(row))
    f.close()
    mtw.finalize()


def cpmd_ener_to_tracks(filename, destination, sub=slice(None), clear=True):
    """Convert a cp2k energy file into separate tracks."""
    names = ["step", "fict_kinectic_energy", "temperature", "potential_energy", "classical_energy", "hamiltonian_energy", "ms_displacement"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtypes = [int, float, float, float, float, float, float]
    dtype = numpy.dtype([  (name, t, 1) for name, t in zip(names, dtypes)  ])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    f = file(filename)
    for line in itertools.islice(f, sub.start, sub.stop, sub.step):
        row = tuple(float(word) for word in line.split()[:7])
        mtw.dump_row(row)
    f.close()
    mtw.finalize()


def cp2k_cell_to_tracks(filename, destination, sub=slice(None), clear=True):
    names = ["cell.a.x", "cell.b.x", "cell.c.x", "cell.a.y", "cell.b.y", "cell.c.y", "cell.a.z", "cell.b.z", "cell.c.z", "cell.a", "cell.b", "cell.c", "cell.alpha", "cell.beta", "cell.gamma"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtype = numpy.dtype([("cell", float, (3,3)),("norms", float, 3),("angles", float, 3)])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    cr = CellReader(filename)
    for cell in itertools.islice(cr, sub.start, sub.stop, sub.step):
        norms = numpy.sqrt((cell**2).sum(axis=0))
        alpha = numpy.arccos(numpy.clip(numpy.dot(cell[:,1],cell[:,2])/norms[1]/norms[2], -1,1))
        beta = numpy.arccos(numpy.clip(numpy.dot(cell[:,2],cell[:,0])/norms[2]/norms[0], -1,1))
        gamma = numpy.arccos(numpy.clip(numpy.dot(cell[:,0],cell[:,1])/norms[0]/norms[1], -1,1))
        mtw.dump_row((cell, norms, (alpha, beta, gamma)))
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

    shape = (len(atom_indexes), 6)
    dtype = numpy.dtype([("cor", float, shape)])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)

    f = file(filename)
    counter = 0
    row = 0
    frame = numpy.zeros(shape, float)
    for line in f:
        if row < len(atom_indexes) or atom_indexes[row] == counter:
            words = line.split()[1:]
            for col, word in enumerate(words):
                frame[row,col] = float(word)
            row += 1
        counter += 1
        if counter == num_atoms:
            mtw.dump_row((frame,))
            counter = 0
            row = 0
    f.close()
    mtw.finalize()


def tracks_to_xyz(prefix, destination, symbols, sub=slice(None), file_unit=angstrom, atom_indexes=None, unit_cell_iter=None, groups=None):
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
    dtype = numpy.dtype([("cor", float, (len(atom_indexes), 3))])
    mtr = MultiTracksReader(filenames, dtype, sub=sub)
    for row in mtr:
        coordinates = row["cor"]
        if unit_cell_iter is not None:
            try:
                uc = unit_cell_iter.next()
            except StopIteration:
                raise ValueError("Not enough frames in the unit cell tracks.")
            if groups is None:
                coordinates -= numpy.dot(uc.cell, numpy.floor(numpy.dot(uc.cell_reciproke, coordinates.transpose()))).transpose()
            else:
                for group in groups:
                    center = coordinates[group].mean(axis=0)
                    coordinates[group] -= numpy.dot(uc.cell, numpy.floor(numpy.dot(uc.cell_reciproke, center)))
        xyz_writer.dump("None", coordinates)
    f.close()


def atrj_to_tracks(filename, destination, sub=slice(None), atom_indexes=None, clear=True):
    """Convert an xyz file into separate tracks."""
    atrj_reader = ATRJReader(filename, sub)

    if atom_indexes is None:
        atom_indexes = range(atrj_reader.num_atoms)
    else:
        atom_indexes = list(atom_indexes)

    filenames = []
    fields = []

    for index in atom_indexes:
        for cor in ["x", "y", "z"]:
            filenames.append(os.path.join(destination, "atom.pos.%07i.%s" % (index, cor)))
    fields.append( ("cor", float, (len(atom_indexes),3)) )
    filenames.append(os.path.join(destination, "time"))
    fields.append( ("time", float, 1) )
    filenames.append(os.path.join(destination, "step"))
    fields.append( ("step", int, 1) )
    filenames.append(os.path.join(destination, "total_energy"))
    fields.append( ("tote", float, 1) )

    dtype = numpy.dtype(fields)
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    for frame in atrj_reader:
        mtw.dump_row((
           frame.coordinates[atom_indexes],
           frame.time, frame.step, frame.total_energy
        ))
    mtw.finalize()


