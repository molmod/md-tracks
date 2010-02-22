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


import shlex
from subprocess import Popen, PIPE


__all__ = ["WrapperError", "Wrapper"]


class WrapperError(Exception):
    pass


class Wrapper(object):
    verbose = False

    def __init__(self, name):
        self.name = name

    def __call__(self, *arguments):
        def iter_args():
            for argument in arguments:
                if isinstance(argument, list):
                    for arg in argument:
                        yield arg
                else:
                    yield argument

        args = []
        for arg in iter_args():
            if (isinstance(arg, float) or isinstance(arg, int)) and arg < 0:
                args.append("'(%s)'" % arg)
            else:
                args.append(str(arg))

        command = "%s %s" % (self.name, " ".join(args))
        if self.verbose:
            print command
        p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, error = p.communicate()
        retcode = p.returncode
        if retcode != 0 or self.verbose:
            print "Command output:"
            print output
            print "Command error:"
            print error
            if retcode != 0:
                raise WrapperError("An error occured while executing command (retcode=%i): \n%s" % (retcode, command))
        return output.split("\n")[:-1]

# when scripts are added, this list must be updated
names = [
    "tr-ac", "tr-ac-error", "tr-ac-fft", "tr-angular-momentum",
    "tr-blav", "tr-calc", "tr-corr", "tr-cwt", "tr-derive", "tr-fit-geom",
    "tr-fit-peaks", "tr-fluct", "tr-from-atrj", "tr-from-cp2k-cell",
    "tr-from-cp2k-ener", "tr-from-cp2k-stress", "tr-from-cpmd-ener",
    "tr-from-cpmd-traj", "tr-from-dlpoly-hist", "tr-from-dlpoly-output",
    "tr-from-lammps-dump", "tr-from-gro",
    "tr-from-txt", "tr-from-xyz", "tr-hist", "tr-ic-bend",  "tr-ic-dihed",
    "tr-ic-dist", "tr-ic-dtl", "tr-ic-oop", "tr-ic-psf", "tr-ic-puckering",
    "tr-integrate", "tr-irfft", "tr-length", "tr-mean-std", "tr-msd",
    "tr-msd-fit", "tr-norm", "tr-pca", "tr-pca-geom", "tr-plot", "tr-qh-entropy", "tr-rdf",
    "tr-reduce", "tr-rfft", "tr-select", "tr-select-rings",
    "tr-shortest-distance", "tr-slice", "tr-spectrum", "tr-split-com",
    "tr-to-txt", "tr-to-xyz", "tr-to-xyz-mode",
]

for name in names:
    python_name = name.replace("-", "_")
    exec("%s = Wrapper(name)" % python_name)
    __all__.append(python_name)


