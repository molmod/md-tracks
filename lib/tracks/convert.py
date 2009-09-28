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


from tracks.core import MultiTracksReader, MultiTracksWriter
from molmod.io import XYZReader, ATRJReader, DLPolyHistoryReader, \
    DLPolyOutputReader, LAMMPSDumpReader, GroReader, XYZWriter
from molmod.units import angstrom, femtosecond, deg, amu, picosecond, bar

import os, numpy, itertools


__all__ = [
    "xyz_to_tracks", "cp2k_ener_to_tracks", "cpmd_ener_to_tracks",
    "cp2k_cell_to_tracks", "cp2k_stress_to_tracks", "cpmd_traj_to_tracks",
    "tracks_to_xyz", "atrj_to_tracks", "dlpoly_history_to_tracks",
    "dlpoly_output_to_tracks",
]


def iter_real_lines(f):
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
    mtw.finish()


def cp2k_ener_to_tracks(filename, destination, sub=slice(None), clear=True):
    """Convert a cp2k energy file into separate tracks."""
    names = ["step", "time", "kinetic_energy", "temperature", "potential_energy", "conserved_quantity"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtypes = [int, float, float, float, float, float]
    dtype = numpy.dtype([  (name, t, 1) for name, t in zip(names, dtypes)  ])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    f = file(filename)
    for line in itertools.islice(iter_real_lines(f), sub.start, sub.stop, sub.step):
        row = [float(word) for word in line.split()[:6]]
        row[1] = row[1]*femtosecond
        mtw.dump_row(tuple(row))
    f.close()
    mtw.finish()


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
    mtw.finish()


def cp2k_cell_to_tracks(filename, destination, sub=slice(None), clear=True):
    names = ["step", "time", "cell.a.x", "cell.a.y", "cell.a.z", "cell.b.x", "cell.b.y", "cell.b.z", "cell.c.x", "cell.c.y", "cell.c.z", "volume", "cell.a", "cell.b", "cell.c", "cell.alpha", "cell.beta", "cell.gamma"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtype = numpy.dtype([("step", int),("time", float),("cell", float, (3,3)),("volume", float),("norms", float, 3),("angles", float, 3)])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    f = file(filename)
    for line in itertools.islice(iter_real_lines(f), sub.start, sub.stop, sub.step):
        values = [float(word) for word in line.split()[:12]]
        row = [int(values[0]),values[1]*femtosecond]
        cell = numpy.array(values[2:11]).reshape(3,3).transpose()*angstrom
        row.append(cell)
        row.append(values[11] * angstrom**3)
        norms = numpy.sqrt((cell**2).sum(axis=0))
        row.append(norms)
        alpha = numpy.arccos(numpy.clip(numpy.dot(cell[:,1],cell[:,2])/norms[1]/norms[2], -1,1))
        beta = numpy.arccos(numpy.clip(numpy.dot(cell[:,2],cell[:,0])/norms[2]/norms[0], -1,1))
        gamma = numpy.arccos(numpy.clip(numpy.dot(cell[:,0],cell[:,1])/norms[0]/norms[1], -1,1))
        row.append(numpy.array([alpha,beta,gamma]))
        mtw.dump_row(tuple(row))
    f.close()
    mtw.finish()


def cp2k_stress_to_tracks(filename, destination, sub=slice(None), clear=True):
    names = ["step", "time", "stress.xx", "stress.xy", "stress.xz", "stress.yx", "stress.yy", "stress.yz", "stress.zx", "stress.zy", "stress.zz", "pressure"]
    filenames = list(os.path.join(destination, name) for name in names)
    dtype = numpy.dtype([("step", int),("time", float),("stress", float, (3,3)),("pressure", float)])
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    f = file(filename)
    for line in itertools.islice(iter_real_lines(f), sub.start, sub.stop, sub.step):
        values = [float(word) for word in line.split()[:11]]
        row = [int(values[0]),values[1]*femtosecond]
        cell = numpy.array(values[2:11]).reshape(3,3).transpose()*bar
        row.append(cell)
        row.append((cell[0,0]+cell[1,1]+cell[2,2])/3)
        mtw.dump_row(tuple(row))
    f.close()
    mtw.finish()


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
    mtw.finish()


def tracks_to_xyz(prefix, destination, symbols, sub=slice(None), file_unit=angstrom, atom_indexes=None, unit_cell_iter=None, groups=None):
    """Converts a set of tracks into an xyz file."""
    if atom_indexes is None:
        atom_indexes = range(len(symbols))
    else:
        atom_indexes = list(atom_indexes)
        if groups is not None:
            # reduce the groups to the selected atoms and use the index of the
            # reduced set.
            reverse_indexes = dict((atom_index, counter) for counter, atom_index in enumerate(atom_indexes))
            new_groups = []
            for group in groups:
                new_group = []
                for atom_index in group:
                    new_index = reverse_indexes.get(atom_index)
                    if new_index is not None:
                        new_group.append(new_index)
                if len(new_group) > 0:
                    new_groups.append(new_group)
            groups = new_groups
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
                coordinates -= numpy.dot(uc.matrix, numpy.floor(numpy.dot(uc.reciprocal, coordinates.transpose()))).transpose()
            else:
                for group in groups:
                    center = coordinates[group].mean(axis=0)
                    coordinates[group] -= numpy.dot(uc.matrix, numpy.floor(numpy.dot(uc.reciprocal, center)))
        xyz_writer.dump("None", coordinates)
    f.close()


def atrj_to_tracks(filename, destination, sub=slice(None), atom_indexes=None, clear=True):
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
    mtw.finish()


def dlpoly_history_to_tracks(
    filename, destination, sub=slice(None), atom_indexes=None, clear=True,
    pos_unit=angstrom, vel_unit=angstrom/picosecond, frc_unit=amu*angstrom/picosecond**2, time_unit=picosecond,
    mass_unit=amu
):
    hist_reader = DLPolyHistoryReader(filename, sub, pos_unit, vel_unit, frc_unit, time_unit, mass_unit)

    if atom_indexes is None:
        atom_indexes = range(hist_reader.num_atoms)
    else:
        atom_indexes = list(atom_indexes)

    filenames = []
    fields = []

    filenames.append(os.path.join(destination, "step"))
    fields.append( ("step", int, 1) )
    filenames.append(os.path.join(destination, "time"))
    fields.append( ("time", float, 1) )
    for vec in "abc":
        for cor in "xyz":
            filenames.append(os.path.join(destination, "cell.%s.%s" % (vec, cor)))
    fields.append( ("cell", float, (3,3)) )
    for vec in "abc":
        filenames.append(os.path.join(destination, "cell.%s" % (vec)))
    fields.append( ("norms", float, 3) )
    for angle in "alpha", "beta", "gamma":
        filenames.append(os.path.join(destination, "cell.%s" % (angle)))
    fields.append( ("angles", float, 3) )
    for index in atom_indexes:
        for cor in "xyz":
            filenames.append(os.path.join(destination, "atom.pos.%07i.%s" % (index, cor)))
    fields.append( ("pos", float, (len(atom_indexes),3)) )
    if hist_reader.keytrj > 0:
        for index in atom_indexes:
            for cor in "xyz":
                filenames.append(os.path.join(destination, "atom.vel.%07i.%s" % (index, cor)))
        fields.append( ("vel", float, (len(atom_indexes),3)) )
    if hist_reader.keytrj > 1:
        for index in atom_indexes:
            for cor in "xyz":
                filenames.append(os.path.join(destination, "atom.frc.%07i.%s" % (index, cor)))
        fields.append( ("frc", float, (len(atom_indexes),3)) )

    dtype = numpy.dtype(fields)
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    for frame in hist_reader:
        cell = frame["cell"]
        norms = numpy.sqrt((cell**2).sum(axis=0))
        frame["norms"] = norms
        alpha = numpy.arccos(numpy.clip(numpy.dot(cell[:,1],cell[:,2])/norms[1]/norms[2], -1,1))
        beta = numpy.arccos(numpy.clip(numpy.dot(cell[:,2],cell[:,0])/norms[2]/norms[0], -1,1))
        gamma = numpy.arccos(numpy.clip(numpy.dot(cell[:,0],cell[:,1])/norms[0]/norms[1], -1,1))
        frame["angles"] = [alpha, beta, gamma]
        frame["pos"] = frame["pos"][atom_indexes]
        if hist_reader.keytrj > 0:
            frame["vel"] = frame["vel"][atom_indexes]
        if hist_reader.keytrj > 0:
            frame["frc"] = frame["frc"][atom_indexes]
        mtw.dump_row(tuple(frame[name] for name, type, shape in fields))
    mtw.finish()


def dlpoly_output_to_tracks(
    filename, destination, sub=slice(None), clear=True, skip_equi_period=True,
    pos_unit=angstrom, time_unit=picosecond, angle_unit=deg, e_unit=amu/(angstrom/picosecond)**2
):
    output_reader = DLPolyOutputReader(filename, sub, skip_equi_period, pos_unit, time_unit, angle_unit, e_unit)

    filenames = [
        "step", "conserved_quantity", "temperature", "potential_energy",
        "vanderwaals_energy", "coulomb_energy", "bond_energy", "bending_energy",
        "torsion_energy", "tethering_energy", "time", "enthalpy",
        "rotational_temperature", "virial", "vanderwaals_virial",
        "coulomb_virial", "bond_viral", "bending_virial", "constraint_virial",
        "tethering_virial", "cputime", "volume", "shell_temperature",
        "shell_energy", "shell_virial", "cell.alpha", "cell.beta", "cell.gamma",
        "pmf_virial", "pressure",
    ]
    filenames = [os.path.join(destination, filename) for filename in filenames]
    fields = [("step", int)] + [("foo%i" % i, float) for i in xrange(29)]

    dtype = numpy.dtype(fields)
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    for row in output_reader:
        mtw.dump_row(tuple(row))
    mtw.finish()


def lammps_dump_to_tracks(filename, destination, meta, sub=slice(None), clear=True):

    units = []
    for unit, name, isvector in meta:
        if isvector:
            units.extend([unit, unit, unit])
        else:
            units.append(unit)

    dump_reader = LAMMPSDumpReader(filename, units, sub)
    num_atoms = dump_reader.num_atoms

    filenames = [os.path.join(destination, "step")]
    fields = [("step", int)]

    for unit, name, isvector in meta:
        if isvector:
            for i in xrange(num_atoms):
                filenames.append(os.path.join(destination, "atom.%s.%07i.x" % (name, i)))
            fields.append(("atom.%s.x" % name, float, num_atoms))
            for i in xrange(num_atoms):
                filenames.append(os.path.join(destination, "atom.%s.%07i.y" % (name, i)))
            fields.append(("atom.%s.y" % name, float, num_atoms))
            for i in xrange(num_atoms):
                filenames.append(os.path.join(destination, "atom.%s.%07i.z" % (name, i)))
            fields.append(("atom.%s.z" % name, float, num_atoms))
        else:
            for i in xrange(num_atoms):
                filenames.append(os.path.join(destination, "atom.%s.%07i" % (name, i)))
            fields.append(("atom.%s" % name, float, num_atoms))


    dtype = numpy.dtype(fields)
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    for frame in dump_reader:
        mtw.dump_row(tuple(frame))
    mtw.finish()


def gro_to_tracks(filename, destination, sub=slice(None), clear=True):
    gro_reader = GroReader(filename, sub)
    num_atoms = gro_reader.num_atoms

    names = ["time"]
    fields = [("time", numpy.float32)]
    names.extend(sum((
        ["atom.pos.%07i.x" % index, "atom.pos.%07i.y" % index,
         "atom.pos.%07i.z" % index]
        for index in xrange(num_atoms)
    ), []))
    fields.append(("pos", numpy.float32, (num_atoms,3)))
    names.extend(sum((
        ["atom.vel.%07i.x" % index, "atom.vel.%07i.y" % index,
         "atom.vel.%07i.z" % index]
        for index in xrange(num_atoms)
    ), []))
    fields.append(("vel", numpy.float32, (num_atoms,3)))
    names.extend([
        "cell.a.x", "cell.b.x", "cell.c.x",
        "cell.a.y", "cell.b.y", "cell.c.y",
        "cell.a.z", "cell.b.z", "cell.c.z",
    ])
    fields.append(("cell", numpy.float32, (3,3)))

    dtype = numpy.dtype(fields)
    filenames = [os.path.join(destination, name) for name in names]
    mtw = MultiTracksWriter(filenames, dtype, clear=clear)
    for time, pos, vel, cell in gro_reader:
        mtw.dump_row((time, pos, vel, cell))
    mtw.finish()

