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


import os, sys, unittest, shutil


__all__ = [
    "orig_dir", "scripts_dir", "lib_dir", "tmp_dir", "input_dir", "output_dir",
    "BaseTestCase",
]


orig_dir = os.getcwd()
scripts_dir = os.path.join(os.path.dirname(os.getcwd()), "scripts")
lib_dir = os.path.join(os.path.dirname(os.getcwd()), "lib")
tmp_dir = os.path.join(os.getcwd(), "tmp")
input_dir = os.path.join(os.getcwd(), "input")
output_dir = os.path.join(os.getcwd(), "output")
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

sys.path.insert(0, lib_dir)


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)
        os.chdir(tmp_dir)

    def tearDown(self):
        os.chdir(orig_dir)
        shutil.rmtree(tmp_dir)
