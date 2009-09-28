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


from molmod.transformations import superpose

import numpy


__all__ = ["fit_geometry"]


def fit_geometry(ref_coordinates, mtr, mtw, do_geom=False, do_transform=False, weights=None):
    """Fits a trajectorie of cartesian coordinates to a reference geometry

       This is the well-known rmsd fit typically performed on a trajectory after
       a long solute/solvent simulation (on the solute). It is based on the
       Kabsch algorithm:
       http://en.wikipedia.org/wiki/Kabsch_algorithm
       http://dx.doi.org/10.1107/S0567739476001873

       This implementation has a few features that are sorely missing in
       comparable software tools. One can optionally use the atomic masses as
       weights, which is a better approximation when one is interested in
       removing globale linear and angular momentum.

       Arguments:
         ref_coordinates  --  the reference coordinates to fit to (Nx3 array)
         mtr  --  A MultiTracksReader that iterates of the frames of the
                  trajectory.
         mtw  --  A MultiTracksWriter that writes out the rmsd, and optionally
                  the rotated and translated geometry (do_geom=True), and the
                  rotation matrix and translation vector. (do_transform=True)
         do_geom  --  When True, the rotated and translated geometry is also
                      written out
         do_transform  --  When True, the rotation matrix and translation vector
                           are also written
         weights  --  The weights to be used in the Kabsch algorithm
    """
    for frame in mtr:
        transform = superpose(ref_coordinates, frame[0], weights)
        new_coordinates = numpy.dot(frame[0], transform.r.transpose()) + transform.t
        rmsd = ((new_coordinates - ref_coordinates)**2).mean()
        row = [rmsd]
        if do_geom:
            row.append(new_coordinates)
        if do_transform:
            row.append(transform.t)
            row.append(transform.r)
        mtw.dump_row(tuple(row))
    mtw.finish()

