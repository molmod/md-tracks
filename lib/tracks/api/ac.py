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


from tracks.api.spectrum import SpectrumProcessor

import numpy


__all__ = [
    "fit_cor_time", "compute_ac_fft", "cor_time", "mean_error_ac",
    "compute_blav", "mean_error_blav"
]


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


def integrate_cor_time(time_step, ac):
    """Integrate the correlation time from a normalized autocorrelation function.

       The integration is carried out over the first positive part of the auto-
       correlation function.

       Arguments:
         time_step  --  The time step in atomic units.
         ac  --  An array with the normalized autocorrelation function.
    """
    tmp = (ac<0).nonzero()[0]
    if len(tmp) == 0:
        return time_step
    end = tmp[0]
    if end <= 1:
        return time_step
    return time_step*ac[:end].sum()


def compute_ac_fft(inputs, num_blocks=10, zero_mean=False):
    """Compute the autocorrelation of a given set of signals

       Argument:
         inputs  --  Arrays of equal length with time dependent information.

       Optional arguments:
         num_blocks  --  The number of blocks used in the computation of the
                         spectrum. [default=10]
         zero_mean  --  When True, the average is not subtracted from the
                        signal.
    """
    sp = SpectrumProcessor(1.0, num_blocks)
    for inp in inputs:
        if not zero_mean:
            inp = inp - inp.mean()
        sp.process(inp)
    amp = sp.get_results()[2]
    return numpy.fft.irfft(amp)


def cor_time(time_step, inputs, num_blocks=10, zero_mean=False):
    """Determine the correlation time for a list of inputs.

       The power spectrum of the inputs is computed and the first part of the
       inverse transform is used to fit a correlation time.

       Arguments:
         time_step --  The time step in atomic units.
         inputs  --  Arrays of equal length with time dependent information.

       Optional arguments
         num_blocks  --  The number of blocks used in the computation of the
                         spectrum. [default=10]
         zero_mean  --  When True, the average is not subtracted from the
                        signal.
    """
    ac = compute_ac_fft(inputs, num_blocks, zero_mean)
    ac /= ac[0]
    return integrate_cor_time(time_step, ac)


def mean_error_ac(signal, std=None, num_blocks=10):
    """Compute the mean and the error on the mean.

       The error on the mean takes into account the time correlation in the
       signal. The autocorrelation function is computed with the fft algorithm.
       See Appendix D of the book "Understanding molecular simulation" by Daan
       Frenkel and Berend Smit for more information.

       Argument:
         signal  --  A time dependent function.

       Optional argumnts:
         std  --  The standard deviation on the signal, computed as signal.std()
                  is not given
         num_blocks  --  The number of blocks used in the fast Fourrier
                         transform.
    """
    tau = cor_time(1.0, [signal], num_blocks)
    mean = signal.mean()
    if std is None:
        std = signal.std()
    return mean, std*numpy.sqrt(2*tau/len(signal))


def compute_blav(time_step, signal, min_blocks=100):
    """Analyze the signal with the block average method.

       Arguments:
         time_step  --  the time step of the second argument.
         signal  --  A time dependent function.

       Optional arguments:
         min_blocks  --  The minimum number of blocks to be considered.

       Returns:
         mean  --  the signal mean
         einf  --  the fitted error
         cinf  --  the fitted correlation time
         sinf  --  the fitted statistical inefficiency
         be  --  the convergence of the fitted error
         bc  --  the convergence of the fitted correlation time
         bs  --  the convergence of the fitted statistical inefficiency
         x  --  the block sizes
         e  --  the errors for each block size
         c  --  the correlation times for each block size
         s  --  the statistical inefficiencies for each block size
         l  --  the length of the part used for the fit
    """
    x = [] # block sizes
    e = [] # error on the mean
    c = [] # correlation time
    s = [] # statistical inefficiency

    for block_size in xrange(1, len(signal)/min_blocks):
        n_blocks = len(signal)/block_size
        total_size = n_blocks * block_size
        averages = numpy.array([
            signal[i*block_size:(i+1)*block_size].mean() for i in xrange(n_blocks)
        ], float)

        x.append(block_size*time_step)
        e.append(averages.std()/numpy.sqrt(n_blocks))
        c.append(0.5*averages.var()/signal[:total_size].var()*block_size*time_step)
        s.append(averages.var()/signal[:total_size].var()*block_size)

    x = numpy.array(x)
    e = numpy.array(e)
    c = numpy.array(c)
    s = numpy.array(s)

    l = len(e)*2/3
    if l == 0:
        raise ValueError("Too few blocks to do a proper estimate of the error.")
    select = numpy.arange(len(x)-l, len(x))
    einf, be = numpy.linalg.lstsq(
        numpy.array([numpy.ones(len(select), float), 1/x[select]]).transpose(),
        e[select],
    )[0]

    select = select[numpy.isfinite(c[select])]
    if len(select) > 3:
        (cinf, bc), resids, rank, svals = numpy.linalg.lstsq(
            numpy.array([numpy.ones(len(select), float), 1/x[select]]).transpose(),
            c[select],
        )
        (sinf, bs), resids, rank, svals = numpy.linalg.lstsq(
            numpy.array([numpy.ones(len(select), float), 1/x[select]]).transpose(),
            s[select],
        )
    else:
        cinf = numpy.inf
        sinf = numpy.inf
        bc = 0
        bs = 0

    return signal.mean(), einf, cinf, sinf, be, bc, bs, x, e, c, s, l


def mean_error_blav(signal, min_blocks=100):
    """Compute the mean and the error on the mean with the block average method.

       Argument:
         signal  --  A time dependent function.

       Optional argument:
         min_blocks  --  The minimum number of blocks to be considered.
    """
    return compute_blav(1.0, signal, min_blocks=min_blocks)[:2]


