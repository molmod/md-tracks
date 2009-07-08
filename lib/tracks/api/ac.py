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


import numpy


__all__ = ["fit_cor_time"]

def fit_cor_time(time_step, ac):
    """Fit the correlation time from a normalized autocorrelation function

       The fit is performed based on a simple single exponential decay model.
       Only the beginning of the ac is used for the fit. As soon as the ac drops
       below 0.4, the rest of the data is discarded.

       Arguments:
         time_step  --  The time step in atomic units
         ac  --  An array with the normalized autocorrelation function
    """
    tmp = (ac<0.4).nonzero()[0]
    if len(tmp) == 0:
        return numpy.nan
    end = tmp[0]
    time = numpy.arange(end)*time_step
    ac_log = numpy.log(ac[:end])
    time_constant = -1.0/(numpy.dot(ac_log[:end], time[:end])/numpy.dot(time[:end], time[:end]))
    return time_constant

