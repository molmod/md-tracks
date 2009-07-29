#!/usr/bin/python
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

from tracks.api import *

import unittest, numpy


__all__ = ["ACTestCase", "PCATestCase"]


class ACTestCase(BaseTestCase):
    def test_mean_error(self):
        length = 1000
        signal = numpy.random.normal(0,0.3,length)
        signal += numpy.cos(numpy.arange(length, dtype=float)/length*10*numpy.pi)
        signal += 0.3
        mean, error = mean_error_ac(signal)
        self.assert_(abs(mean - 0.3) < 0.1)
        self.assert_(error < 0.15)
        mean, error = mean_error_blav(signal)
        self.assert_(abs(mean - 0.3) < 0.1)
        self.assert_(error < 0.15)


class PCATestCase(BaseTestCase):
    def check_sanity(self, overlap_fn, num):
        size = 10
        trials = 10
        for trial in xrange(trials):
            # construct 'num' random positive definite matrices
            covs = []
            for i in xrange(num):
                A = numpy.random.normal(0, 1, (size,size))
                covs.append(numpy.dot(A.transpose(), A))

            # test the scale invariance of the overlap function
            scale = numpy.random.uniform(1,2)
            covs_mod = [cov*scale for cov in covs]
            self.assertAlmostEqual(overlap_fn(covs), overlap_fn(covs_mod))

            # test the invariance of the overlap function under a unitary
            # tranformation
            tmp = numpy.random.normal(0, 1, (size,size))
            U,S,Vt = numpy.linalg.svd(tmp)
            del S
            del Vt
            covs_mod = [numpy.dot(U.transpose(), numpy.dot(cov, U)) for cov in covs]
            self.assertAlmostEqual(overlap_fn(covs), overlap_fn(covs_mod))

    def test_cov_overlap(self):
        def my_overlap(covs):
            return cov_overlap(*covs)
        self.check_sanity(my_overlap, 2)

    def test_cov_overlap_multi(self):
        self.check_sanity(cov_overlap_multi, 3)


