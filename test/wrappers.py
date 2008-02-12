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
# Contact information:
#
# Supervisors
#
# Prof. Dr. Michel Waroquier and Prof. Dr. Ir. Veronique Van Speybroeck
#
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 59
# Fax: +32 9 264 65 60
# Email: Michel.Waroquier@UGent.be
# Email: Veronique.VanSpeybroeck@UGent.be
#
# Author
#
# Ir. Toon Verstraelen
# Center for Molecular Modeling
# Ghent University
# Proeftuinstraat 86, B-9000 GENT - BELGIUM
# Tel: +32 9 264 65 56
# Email: Toon.Verstraelen@UGent.be
#
# --


from common import *

from tracks.wrappers import tr_read
from tracks.core import dump_track
from tracks.log import log

import unittest, numpy, os, glob


log.verbose = False


class WrapperTestCase(BaseTestCase):
    def test_tr_read(self):
        values = numpy.arange(0, 10, 0.1, float)
        dump_track(os.path.join(tmp_dir, "tmp"), values)
        output = tr_read(os.path.join(tmp_dir, "tmp"))
        check_values = numpy.array([float(word) for word in output], float)
        self.assertArraysAlmostEqual(values, check_values, 1e-10)

    def test_names(self):
        from tracks.wrappers import names
        names = set(names)
        for name in names:
            self.assert_(os.path.isfile(os.path.join(scripts_dir, name)), "%s missing in scripts_dir" % name)
        for filename in glob.glob(os.path.join(scripts_dir, "*")):
            self.assert_((os.path.basename(filename) in names), "%s missing in names list" % filename)

