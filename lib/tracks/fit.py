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


import numpy


__all__ = ["Model", "GaussianModel", "PeakModel", "FitCostFunction"]


class Model(object):
    def get_parameters(self):
        raise NotImplementedError

    def set_parameters(self, parameters):
        raise NotImplementedError

    def __call__(self, f):
        raise NotImplementedError

    def gradient(self, f):
        raise NotImplementedError


class GaussianModel(Model):
    def __init__(self, parameters):
        assert parameters.shape == (3,)
        self.parameters = parameters

    def get_parameters(self):
        return tuple(self.parameters)

    def set_parameters(self, parameters):
        assert parameters.shape == (3,)
        self.parameters = parameters.copy()

    def get_labels(self):
        return 'lc', 'sigma', 'A'

    def __call__(self, l):
        lc, sigma, A = self.get_parameters()
        return A*numpy.exp(-((l - lc)/sigma)**2)

    def gradient(self, l):
        lc, sigma, A = self.get_parameters()
        x = (l - lc)/sigma
        g = numpy.exp(-x**2)
        return numpy.array([
            A*g*2*x/sigma,
            A*g*2*x**2/sigma,
            g,
        ], float)

    def hessian(self, l):
        lc, sigma, A = self.get_parameters()
        x = (l - lc)/sigma
        g = numpy.exp(-x**2)
        result = numpy.zeros((3,3,len(l)),float)
        result[0,0] = A*g*2/sigma**2*(2*x**2-1)
        result[0,1] = A*g*4*x/sigma**2*(x**2-1)#A*g*2*(x/sigma)**2*(2*x-3)
        result[1,0] = result[0,1]
        result[0,2] = g*2*x/sigma
        result[2,0] = result[0,2]
        result[1,1] = A*g*2*(x/sigma)**2*(2*x**2-3)
        result[1,2] = g*2*x**2/sigma
        result[2,1] = result[1,2]
        result[2,2] = 0
        return result


class PeakModel(Model):
    def __init__(self, parameters=None):
        if parameters is None:
            self.b = 0
            self.gaussians = [GaussianModel()]
        else:
            self.b = parameters[0]
            num_gaussians = (len(parameters)-1)/3
            self.gaussians = [GaussianModel(parameters[index*3+1:index*3+4]) for index in xrange(num_gaussians)]

    def get_parameters(self):
        result = (self.b,)
        for gaussian in self.gaussians:
            result += gaussian.get_parameters()
        return result

    def set_parameters(self, parameters):
        self.b = parameters[0]
        for index, gaussian in enumerate(self.gaussians):
            gaussian.set_parameters(parameters[1+index*3:4+index*3])

    def get_labels(self):
        return ('b',) + sum((gaussian.get_labels() for gaussian in self.gaussians), ())

    def __call__(self, l):
        return self.b+sum(gaussian(l) for gaussian in self.gaussians)

    def gradient(self, l):
        return numpy.array([
            numpy.ones(len(l), float)
        ] + sum(
            (list(gaussian.gradient(l)) for gaussian in self.gaussians)
        , []))

    def hessian(self, l):
        num_parameters = len(self.get_parameters())
        result = numpy.zeros((num_parameters,num_parameters,len(l)), float)
        for index, gaussian in enumerate(self.gaussians):
            result[index*3+1:index*3+4,index*3+1:index*3+4] = gaussian.hessian(l)
        return result


class FitCostFunction(object):
    def __init__(self, wavenumbers, amplitudes, peak_model):
        self.wavenumbers = wavenumbers
        self.amplitudes = amplitudes
        self.peak_model = peak_model

    def __call__(self, parameters):
        self.peak_model.set_parameters(parameters)
        return ((self.peak_model(self.wavenumbers)-self.amplitudes)**2).mean()

    def gradient(self, parameters):
        self.peak_model.set_parameters(parameters)
        return 2*numpy.dot(self.peak_model.gradient(self.wavenumbers), (self.peak_model(self.wavenumbers)-self.amplitudes))/len(self.wavenumbers)

    def hessian(self, parameters):
        errors = (self.peak_model(self.wavenumbers) - self.amplitudes)
        mg = self.peak_model.gradient(self.wavenumbers)
        mh = self.peak_model.hessian(self.wavenumbers)
        result = numpy.zeros((len(parameters),len(parameters)), float)
        for index1 in xrange(len(parameters)):
            for index2 in xrange(len(parameters)):
                result[index1,index2] = numpy.dot(mg[index1],mg[index2]) + numpy.dot(errors,mh[index1,index2])
        return 2*result/len(self.wavenumbers)

