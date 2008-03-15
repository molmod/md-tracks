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


from common import *

from tracks.core import dump_track
from tracks.parse import *
from tracks.log import log

import unittest, numpy


log.verbose = False


__all__ = ["ParseTestCase"]


class ParseTestCase(BaseTestCase):
    def test_step(self):
        dump_track("test", numpy.arange(50))
        self.assertEqual(parse_x_step("test"), 1)

    def test_duration(self):
        dump_track("test", numpy.arange(50))
        self.assertEqual(parse_x_duration("test"), 49)

    def test_length(self):
        dump_track("test", numpy.arange(50))
        self.assertEqual(parse_x_length("test"), 50)

