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


from molmod.constants import lightspeed

import numpy


__all__ = ["compute_spectrum", "SpectrumProcessor"]


def compute_spectrum(time_step, inputs, num_blocks=10):
    """A simple interface to the SpectrumProccessor.

       Arguments:
         time_step  --  The time step [au] of the inputs that will be given
                        to the process method.
         inputs  --  A list of one-dimensional arrays with function values taken
                     at equi-distant time steps.
         num_blocks  --  The number of blocks in which the inputs will be
                         divided. A Fourier transform of each block is computed
                         and the final spectrum is the average over all blocks.
                         [default=10] A larger number of blocks implies a lower
                         resolution on the frequency scale, but yields an
                         improved statistical accuracy of the amplitudes.

       Returns a tuple with four elements:
         frequency_res -- The resolution of the spectrum on the frequency
                          scale.
         wavenumber_res  --  The resulution of the spectrum on the
                             wavenumber scale.
         amplitudes  --  An array of amplitudes representing the spectrum.
         amplitudes_err  --  The statistical error on the amplitudes.
    """
    sp = SpectrumProcessor(time_step, num_blocks)
    for i in inputs:
        sp.process(i)
    return sp.get_results()


class SpectrumProcessor(object):
    """Implements a stepwise implementation for the computation of spectra.

       Each input from a multidimensional trajectory is processed separately.
       One does not have to load the entire trajectory in memory at the same
       time.
    """

    def __init__(self, time_step, num_blocks=10):
        """Initialize a spectrum SpectrumProcessor instance.

           Arguments:
             time_step  --  The time step [au] of the inputs that will be given
                            to the process method.
             num_blocks  --  The number of blocks in which the inputs will be
                             divided. A Fourier transform of each block is computed
                             and the final spectrum is the average over all blocks.
                             [default=10] A larger number of blocks implies a lower
                             resolution on the frequency scale, but yields an
                             improved statistical accuracy of the amplitudes.
        """
        self._time_step = time_step
        self._num_blocks = num_blocks
        self._sum = 0
        self._sum_sq = 0
        self._count = 0

        self._input_length = None
        self._block_size = None

    block_size = property(lambda self: self._block_size)

    def process(self, fn):
        """Add an input to the spectrum.

           Arguments:
             fn  --  A one-dimensional array with function values taken at equi-
                     distant time steps. (See the time_step argument of the
                     __init__ method.)
        """
        if len(fn) < 2*self._num_blocks:
            raise ValueError("The length of the input for the spectrum must at least be two times the block size.")
        if self._input_length is None:
            self._input_length = len(fn)
            self._block_size = len(fn)/self._num_blocks
        elif self._input_length != len(fn):
            raise ValueError("All inputs must have the same length.")

        for index in xrange(self._num_blocks):
            fn_block = fn[index*self._block_size:(index+1)*self._block_size]
            amp_block = abs(numpy.fft.rfft(fn_block))**2
            self._sum += amp_block
            self._sum_sq += amp_block**2
            self._count += 1

    def get_results(self):
        """Compute the resulting spectrum.

           Returns a tuple with four elements:
             frequency_res -- The resolution of the spectrum on the frequency
                              scale.
             wavenumber_res  --  The resulution of the spectrum on the
                                 wavenumber scale.
             amplitudes  --  An array of amplitudes representing the spectrum.
             amplitudes_err  --  The statistical error on the amplitudes.
        """
        if self._count == 0:
            raise RuntimeError("There are no results yet. Call the process method first.")

        amplitudes = self._sum/(self._count*self._block_size)
        amplitudes_sq = self._sum_sq/(self._count*self._block_size)
        amplitudes_err = numpy.sqrt((amplitudes_sq - amplitudes**2))/self._count # error on the mean

        duration = self._time_step*self._block_size
        frequency_res = 1.0/duration
        wavenumber_res = 1.0/duration/lightspeed

        return frequency_res, wavenumber_res, amplitudes, amplitudes_err


