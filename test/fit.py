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


from tracks.fit import PeakModel, FitCostFunction

import unittest, numpy


__all__ = ["FitTestCase"]


class FitTestCase(unittest.TestCase):
    def setUp(self):
        self.initial_parameters = numpy.array([0.1,560,50,1.0,610,30,0.3], float)
        self.peak_model = PeakModel(self.initial_parameters)
        self.wavenumbers = numpy.arange(400,700,3.37,dtype=float)
        self.target_data = numpy.random.normal(0,1,self.wavenumbers.shape)
        self.cost_function = FitCostFunction(self.wavenumbers, self.target_data, self.peak_model)

    def test_model_gradient(self):
        self.peak_model.set_parameters(self.initial_parameters)
        epsilon = 1e-5
        ag = self.peak_model.gradient(self.wavenumbers) # analytical gradient
        for index in xrange(len(self.initial_parameters)):
            tmp = self.initial_parameters.copy()
            tmp[index] -= 0.5*epsilon
            self.peak_model.set_parameters(tmp)
            m = self.peak_model(self.wavenumbers)
            tmp = self.initial_parameters.copy()
            tmp[index] += 0.5*epsilon
            self.peak_model.set_parameters(tmp)
            p = self.peak_model(self.wavenumbers)
            #print "target = (p-m)/epsilon"
            target = (p-m)/epsilon
            #print target
            #print "outcome = ag[index]"
            outcome = ag[index]
            #print outcome
            #print "ratio = outcome/target"
            #print outcome/target
            error = abs((p-m)/epsilon - ag[index]).max()
            oom = abs(ag[index]).max()
            #print error,oom
            self.assert_(error/oom < 1e-2, "model gradient incorrect, parameter %i" % index)

    def test_model_hessian(self):
        self.peak_model.set_parameters(self.initial_parameters)
        epsilon = 1e-5
        ah = self.peak_model.hessian(self.wavenumbers) # analytical hessian
        for index1 in xrange(len(self.initial_parameters)):
            tmp = self.initial_parameters.copy()
            tmp[index1] -= 0.5*epsilon
            self.peak_model.set_parameters(tmp)
            m = self.peak_model.gradient(self.wavenumbers)
            tmp = self.initial_parameters.copy()
            tmp[index1] += 0.5*epsilon
            self.peak_model.set_parameters(tmp)
            p = self.peak_model.gradient(self.wavenumbers)
            for index2 in xrange(len(self.initial_parameters)):
                #print "index1 index2 =", index1, index2
                #print "target = (p[index2]-m[index2])/epsilon"
                target = (p[index2]-m[index2])/epsilon
                #print target
                #print "outcome = ah[index1,index2]"
                outcome = ah[index1,index2]
                #print outcome
                #print "target/outcome"
                #print target/outcome
                error = abs(outcome - target).max()
                oom = abs(target).max()
                #print error,oom
                if oom == 0:
                    self.assert_(error < epsilon, "model hessian incorrect, parameters (%i,%i)" % (index1,index2))
                else:
                    self.assert_(error/oom < 1e-2, "model hessian incorrect, parameters (%i,%i)" % (index1,index2))

    def test_cost_gradient(self):
        epsilon = 1e-5
        ag = self.cost_function.gradient(self.initial_parameters) # analytical gradient
        for index in xrange(len(self.initial_parameters)):
            tmp = self.initial_parameters.copy()
            tmp[index] -= 0.5*epsilon
            m = self.cost_function(tmp)
            tmp = self.initial_parameters.copy()
            tmp[index] += 0.5*epsilon
            p = self.cost_function(tmp)
            #print "(p-m)/epsilon"
            #print (p-m)/epsilon
            #print "ag[index]"
            #print ag[index]
            #print "ratio"
            #print ag[index]/(p-m)/epsilon
            error = abs((p-m)/epsilon - ag[index]).max()
            oom = abs(ag[index]).max()
            #print error,oom
            self.assert_(error/oom < 1e-2, "model gradient incorrect, parameter %i" % index)

    def test_cost_hessian(self):
        epsilon = 1e-5
        ah = self.cost_function.hessian(self.initial_parameters) # analytical hessian
        for index1 in xrange(len(self.initial_parameters)):
            tmp = self.initial_parameters.copy()
            tmp[index1] -= 0.5*epsilon
            m = self.cost_function.gradient(tmp)
            tmp = self.initial_parameters.copy()
            tmp[index1] += 0.5*epsilon
            p = self.cost_function.gradient(tmp)
            for index2 in xrange(len(self.initial_parameters)):
                #print "index1 index2 =", index1, index2
                #print "target = (p[index2]-m[index2])/epsilon"
                target = (p[index2]-m[index2])/epsilon
                #print target
                #print "outcome = ah[index1,index2]"
                outcome = ah[index1,index2]
                #print outcome
                #print "target/outcome"
                #print target/outcome
                error = abs(outcome - target).max()
                oom = abs(target).max()
                #print error,oom
                if oom == 0:
                    self.assert_(error < epsilon, "cost hessian incorrect, parameters (%i,%i)" % (index1,index2))
                else:
                    self.assert_(error/oom < 1e-2, "cost hessian incorrect, parameters (%i,%i)" % (index1,index2))



