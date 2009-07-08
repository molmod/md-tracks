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


from common import *

from tracks.wrappers import tr_to_txt, tr_rfft
from tracks.core import dump_track
from tracks.log import log

import unittest, numpy, os, glob


__all__ = ["WrapperTestCase"]


log.verbose = False
# make sure the wrappers point to the scripts in the script dir.
tr_to_txt.name = os.path.join(scripts_dir, tr_to_txt.name)
tr_rfft.name = os.path.join(scripts_dir, tr_rfft.name)


class WrapperTestCase(BaseTestCase):
    def test_tr_to_txt(self):
        values = numpy.arange(0, 10, 0.1, float)
        dump_track(os.path.join(tmp_dir, "tmp"), values)
        output = tr_to_txt(-1.0, os.path.join(tmp_dir, "tmp"))
        check_values = numpy.array([-float(word) for word in output], float)
        self.assertArraysAlmostEqual(values, check_values, 1e-10)

    def test_names(self):
        if lib_dir is not None:
            from tracks.wrappers import names
            names = set(names)
            for name in names:
                self.assert_(os.path.isfile(os.path.join(scripts_dir, name)), "%s missing in scripts_dir" % name)
            for filename in glob.glob(os.path.join(scripts_dir, "*")):
                self.assert_((os.path.basename(filename) in names), "%s missing in names list" % filename)

    def test_tr_rfft(self):
        num = 20
        for i in xrange(num):
            dump_track("file.%07i" % i, numpy.random.normal(0,1,100))
        tr_rfft(glob.glob("file.*"))


