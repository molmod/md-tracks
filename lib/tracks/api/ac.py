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


from tracks.api.spectrum import SpectrumProcessor

import numpy


__all__ = ["fit_cor_time", "cor_time", "mean_error"]


def fit_cor_time(time_step, ac):
    """Fit the correlation time from a normalized autocorrelation function.

       The fit is performed based on a simple single exponential decay model.
       Only the beginning of the ac is used for the fit. As soon as the ac drops
       below 1/e, the rest of the data is discarded. If no significant time
       correlation is found, the time_step is returned as correlation time.

       Arguments:
         time_step  --  The time step in atomic units.
         ac  --  An array with the normalized autocorrelation function.
    """
    tmp = (ac<0.3678794).nonzero()[0]
    if len(tmp) == 0:
        return time_step
    end = tmp[0]
    if end <= 1:
        return time_step
    time = numpy.arange(end)*time_step
    ac_log = numpy.log(ac[:end])
    correlation_time = -1.0/(numpy.dot(ac_log[:end], time[:end])/numpy.dot(time[:end], time[:end]))
    return correlation_time


def cor_time(time_step, inputs, num_blocks=10, zero_mean=False):
    """Determine the correlation time for a list of inputs.

       The power spectrum of the inputs is computed and the first part of the
       inverse transform is used to fit a correlation time.

       Arguments:
         time_setp --  The time step in atomic units.
         inputs  --  Arrays of equal length with time dependent information.

       Optional arguments
         num_blocks  --  The number of blocks used in the computation of the
                         spectrum. [default=10]
         zero_mean  --  When True, the average is not subtracted from the
                        signal.
    """
    sp = SpectrumProcessor(time_step, num_blocks)
    for inp in inputs:
        if not zero_mean:
            inp = inp - inp.mean()
        sp.process(inp)
    freq_res, wave_res, amp, amp_err = sp.get_results()
    ac = numpy.fft.irfft(amp)
    ac /= ac[0]
    return fit_cor_time(time_step, ac)


def mean_error(signal, num_blocks=10):
    """Compute the mean and the error on the mean.

       The error on the mean takes into account the time correlation in the
       signal.

       Argument:
         signal  --  A time dependent function

       Optional argumnts:
         num_blocks  --  The number of blocks used in the fast Fourrier
                         transform
    """
    tau = cor_time(1.0, [signal], num_blocks)
    mean = signal.mean()
    std = signal.std()
    return mean, std*numpy.sqrt(tau/len(signal))

