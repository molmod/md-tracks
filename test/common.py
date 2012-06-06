# -*- coding: utf-8 -*-
# MD-Tracks is a trajectory analysis toolkit for molecular dynamics
# and monte carlo simulations.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>, Center
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
#--


import os, sys, glob

orig_dir = os.path.abspath('test')
tmp_dir = os.path.join(orig_dir, 'tmp')
input_dir =  os.path.join(orig_dir, 'input')
output_dir =  os.path.join(orig_dir, 'output')
lib_dir =  os.path.join(orig_dir, '../')
scripts_dir =  os.path.join(orig_dir, '../scripts')

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

import unittest, shutil, numpy


__all__ = [
    "orig_dir", "scripts_dir", "lib_dir", "tmp_dir", "input_dir", "output_dir",
    "BaseTestCase",
]


class BaseTestCase(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)
        os.chdir(tmp_dir)

    def tearDown(self):
        os.chdir(orig_dir)
        shutil.rmtree(tmp_dir)

    def assertArraysEqual(self, a, b):
        self.assertEqual(a.shape, b.shape, "The array shapes do not match.")
        self.assert_((a==b).all(), "The array values do not match.")

    def assertArrayConstant(self, arr, const):
        self.assert_((arr==const).all(), "Some/All array values do not match the constant.")

    def assertArraysAlmostEqual(self, a, b, err_threshold=1e-5, do_abs=False, do_mean=False, verbose=False):
        def log(s):
            if verbose: print s
        if a.shape != b.shape:
            self.fail("Array shapes do not match: %s!=%s" % (a.shape, b.shape))
        if do_abs:
            if do_mean:
                abserr = abs(a-b).mean()
            else:
                abserr = abs(a-b).max()
            log("both")
            log(numpy.hstack([a,b]))
            log("difference")
            log(a-b)
            log("abserr: %s" % abserr)
            if abserr > err_threshold:
                self.fail("The absolute error is too large: %.3e > %.3e" % (abserr, err_threshold))
        else:
            if do_mean:
                relerr = abs(a-b).mean()*2/(abs(a).mean()+abs(b).mean())
            else:
                relerr = abs(a-b).max()*2/(abs(a).max()+abs(b).max())
            log("both")
            log(numpy.hstack([a,b]))
            #log(a)
            #log(b)
            log("difference")
            log(a-b)
            log("relerr: %s" % relerr)
            if relerr > err_threshold:
                self.fail("The relative error is too large: %.3e > %.3e" % (relerr, err_threshold))
            if numpy.isnan(relerr):
                self.fail("The relative error is nan.")
        if numpy.isnan(a).any():
            self.fail("The first argument contains nan's.")
        if numpy.isnan(b).any():
            self.fail("The second argument contains nan's.")

    def assertArrayAlmostConstant(self, arr, const, relerr_threshold):
        error = abs(arr-const).max()
        oom = const
        relerr = error/oom
        self.assert_(relerr <= relerr_threshold, "The relative error is larger than given threshold: %5.3e > %5.3e" % (relerr, relerr_threshold))

    def assertArrayAlmostZero(self, arr, abserr_threshold):
        abserr = abs(arr).max()
        self.assert_(abserr <= abserr_threshold, "The absolute error is larger than given threshold: %5.3e > %5.3e" % (abserr, abserr_threshold))


