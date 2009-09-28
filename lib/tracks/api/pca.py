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


from tracks.core import MultiTracksReader, dump_track, MultiTracksWriter

import numpy


__all__ = [
    "pca_levels", "CovarianceMatrix", "CovarianceBlocks",
    "cov_overlap", "cov_overlap_multi",
    "pca_common_usage", "pca_common_script",
]


def pca_levels(mtr, num_levels, weights=None, correlation=False, reference=None):
    """Perform a full principal component analysis at different block levels

       Argumets:
         mtr  --  A MultiTracksReader that iters of the input at each time step.
         num_levels  -- The number of levels of block sizes to perform the pca
                        on. The first level is all data in one block, the second
                        level has the data divide in two equal blocks, then
                        four equal blocks, then eight equal blocks and so on.

       Optional arguments:
         weights  --  When given, the principal components are computed in
                      weighted coordinates.
         correlation  --  When True, the analysis is performed on the matrix
                          with correlation coefficients. This might be usefull
                          to compare the eigenvalue spectrum with the Wishart
                          distribution.
         reference  --  When given, the reference is assumed to be the vector
                        with the averages of the inputs. Otherwise, the averages
                        are derived from the input.

       Returns:
         cm  --  A covariance matrix object for the entire trajectory. (if the
                 number of frames is uneven, the last frame is dropped.)
         overlap  --  The overlap between the covariance matrix of the first
                      and the last half of the trajectory. (see cov_overlap)

    """
    if weights is not None and correlation:
        raise ValueError("Weighted coordinates have no effect when computing the correlation matrix.")

    block_size = mtr.shortest/(2**(num_levels-1))
    cb_max = CovarianceBlocks(block_size, weights, correlation, reference)
    for data in mtr.iter_buffers():
        data = data["data"]
        cb_max.add_data(data)
    cbs = [cb_max]
    for i in xrange(num_levels-1):
        cbs.insert(0, cbs[0].reduce_blocks())

    return cbs[0].blocks[0], cov_overlap(cbs[1].blocks[0].cov, cbs[1].blocks[1].cov), cbs


class CovarianceMatrix(object):
    """A container for all information related to a covariance matrix."""
    def __init__(self, length, matrix, sum, weights=None, correlation=False, reference=None):
        """Initialize a covariance matrix object

           The arguments to initialize a CovarianceMatrix instance are built up
           by processing a trajectory. Look at CovarianceBlocks for an example.
           They are all stored as attributes. Also some direved properties are
           computed during the initialization. In a second stage, one can
           reprocess that data with init_proj, data_proj and finish_proj to
           compute the principal components and the cosine content.

           Arguments:
             length  --  the length of (a part of) the trajectory used to
                         construct the matrix
             matrix  --  the matrix built up by adding
                         numpy.data(data.transpose(),data) for multiple data
                         arrays belonging to one block
             sum  --  the sum of the data over time

           Optional arguments:
             weights  --  When given, the principal components are computed in
                          weighted coordinates.
             correlation  --  When True, the analysis is performed on the matrix
                              with correlation coefficients. This might be usefull
                              to compare the eigenvalue spectrum with the Wishart
                              distribution.
             reference  --  When given, the reference is assumed to be the vector
                            with the averages of the inputs. Otherwise, the averages
                            are derived from the input.

           Derive attributes:
             cov  --  The actual covariance/correlation matrix
             mean  --  The average of the inputs over time
             evals  --  The eigenvalues of the covariance matrix
             evecs  --  The corresponding eigenvectors
             sigmas  --  The square roots of the eigenvalues.

           Attributes available after projection:
             sqnorms_data  --  squared norms of the principal components
             sqnorms_cos  --  squared norms of the cosines
             dot_data_cos  --  inner product between principal component and
                               cosine
             ccs  --  the cosine contents of each principal component
        """
        # the raw data
        self.length = length
        self.matrix = matrix
        self.sum = sum
        self.weights = weights
        self.correlation = correlation
        self.reference = reference
        # the derived properties
        self.cov = self.matrix / self.length # the actual covariance matrix
        if self.reference is None:
            self.mean = self.sum / self.length
        else:
            self.mean = self.reference
        self.cov -= numpy.outer(self.mean, self.mean)
        if self.correlation:
            diag_sqrt = numpy.sqrt(numpy.diag(self.cov))
            self.cov /= numpy.outer(diag_sqrt, diag_sqrt)
        elif self.weights is not None:
            scale = numpy.sqrt(self.weights)
            self.cov *= numpy.outer(scale, scale)
        # the eigen decomposition
        self.evals, self.evecs = numpy.linalg.eigh(self.cov)
        # largest eigenvalues first
        self.evals = self.evals[::-1]
        self.evecs = self.evecs[:,::-1]
        self.sigmas = numpy.sqrt(abs(self.evals))

    def init_proj(self, output_prefix=None):
        """Setup the projection of the trajectory on the pca eigenmodes.

           When output prefix is given, the principal components are written
           to tracks with the following filenames: ${output_prefix}.${index}.

           After init_proj, call data_proj one or more times with the relevant
           data segments from the trajectory. Finally, call finish_proj.
        """
        N = len(self.evals)
        if output_prefix is not None:
            paths_out = [
                "%s.pc.%07i" % (output_prefix, index)
                for index in xrange(N)
            ]
            dtype = numpy.dtype([("data", float, N)])
            self.proj_mtw = MultiTracksWriter(paths_out, dtype)
        else:
            self.proj_mtw = None
        self.sqnorms_data = numpy.zeros(N, float)
        self.sqnorms_cos = numpy.zeros(N, float)
        self.dot_data_cos = numpy.zeros(N, float)
        self.proj_counter = 0

    def data_proj(self, data):
        """Process data to compute the principal components and the cosine content

           First call init_proj, then call this routine multiple times. Finally
           call finish_proj.
        """
        data = data - self.mean
        if self.correlation:
            data /= numpy.sqrt(numpy.diag(self.cov))
        elif self.weights is not None:
            data *= self.weights
        pcs = numpy.dot(data, self.evecs)
        if self.proj_mtw is not None:
            self.proj_mtw.dump_buffer({"data": pcs})
        t = numpy.arange(self.proj_counter, self.proj_counter+len(data))*(numpy.pi/self.length)
        for i in xrange(data.shape[1]): # iterate ove the columns
            c = numpy.cos((1+i)*t)
            self.sqnorms_data[i] += numpy.dot(pcs[:,i], pcs[:,i])
            self.sqnorms_cos[i] += numpy.dot(c, c)
            self.dot_data_cos[i] += numpy.dot(pcs[:,i], c)
        self.proj_counter += len(data)

    def finish_proj(self):
        """Compute the actual cosine contents.

           Call finish_proj after the last call to data_proj.
        """
        self.ccs = self.dot_data_cos**2/(self.sqnorms_data*self.sqnorms_cos)
        if self.proj_mtw is not None:
            self.proj_mtw.finish()
        del self.proj_mtw
        del self.proj_counter


class CovarianceBlocks(object):
    """A tool to compute covariance matrices efficiently at different block sizes.

       Start with an instance with small blocks sizes, preferentially using the
       following formula:

         blocks_size = total_size / (2**(levels-1))

       Where levels is the number of different block_sizes one wants to
       consider. Then call the process method with data arrays. (Rows correspond
       to different time steps, columns are different coordinates/fields/...)

       A list of covariance matrix objects, self.blocks, is generated gradually
       when more data is provided with the add_data method.

       The reduce method can be used to obtain similar information for large
       block sizes. It takes one argument, num, which is the number of blocks
       will be put toghether to form a new block. A new CovarianceBlocks
       instance is returned as if it was constructed with a larger block_size,
       equal to the orignal block_size times num. When num is larger than the
       number of available blocks, an error is raised.

       Finally one can reprocess all the data with the method project_data to
       project the trajectories on the eigen modes as to compute the principal
       components and the cosine content associated with each mode.
    """

    def __init__(self, block_size, weights=None, correlation=False, reference=None):
        """Intialize a CovarianceBlocks object.

           Arguments:
             block_size  --  The sizes of the data blocks for which covariance
                             matrices are constructed.
             weights  --  The weights associated with the coordinates. When not
                          given, all weights are equal to one.
             correlation  --  Produce matrices with correlation coefficients
                              intead of covariance matrices. (default=False)
             reference  --  A predefined mean for the input signals. If not
                            given the mean is derived from the input data.
        """
        if weights is not None and correlation:
            raise ValueError("Weighted coordinates have no effect when computing the correlation matrix.")
        if block_size <= 0:
            raise ValueError("The second argument, block_size, must be strictly positive.")
        self.block_size = block_size
        self.weights = weights
        self.correlation = correlation
        self.reference = reference

        self.blocks = []
        self._init_matrix()
        self._init_proj()

    def _init_matrix(self):
        """Private methode used by add_data"""
        self._current_length = 0
        self._current_matrix = 0.0
        if self.reference is None:
            self._current_sum = 0.0
        else:
            self._current_sum = None

    def _add_matrix(self, data):
        """Private methode used by add_data"""
        self._current_length += len(data)
        self._current_matrix += numpy.dot(data.transpose(), data)
        if self.reference is None:
            self._current_sum += data.sum(axis=0)

    def _finish_matrix(self):
        """Private methode used by add_data"""
        self.blocks.append(CovarianceMatrix(
            self._current_length, self._current_matrix, self._current_sum,
            self.weights, self.correlation, self.reference
        ))
        self._init_matrix()

    def add_data(self, data):
        """Provide new data to compute covariance matrices.

           Argument:
             data  --  a new array with input data: rows correspond to time
                       frames and columns to coordinates/fields/...

           Make sure the number of columns is always the same for one
           CovarianceBlocks instance.
        """
        pos = 0
        while pos < len(data):
            remaining = len(data) - pos
            todo = self.block_size - self._current_length
            if remaining >= todo:
                self._add_matrix(data[pos:pos+todo])
                pos += todo
                self._finish_matrix()
            else: # todo > remaining
                self._add_matrix(data[pos:])
                break

    def reduce_blocks(self, num=2):
        """Create a new CovarianceBlocks object based on fewer, but larger blocks.

           Argument:
             num  --  The reduction factor for the number of blocks.

           The return CovarianceBlocks object acts as if it was constructed with
           a block size equal to the original block size times the num argument.
           However it is computationally much cheaper to call this method.
        """
        if num > len(self.blocks):
            raise ValueError("Not enough blocks: num=%i > len(self.covs)=%i" % (num, len(self.blocks)))
        result = CovarianceBlocks(self.block_size*num, self.weights, self.correlation, self.reference)
        for i in xrange(len(self.blocks)/num):
            for j in xrange(num):
                result._current_matrix += self.blocks[i*num+j].matrix
                if self.reference is None:
                    result._current_sum += self.blocks[i*num+j].sum
                result._current_length += self.blocks[i*num+j].length
            result._finish_matrix()
        return result

    def _init_proj(self):
        self.proj_counter = -1
        self.proj_done = 0

    def project_data(self, data, output_prefix=None):
        """Process the data in a second run to compute the principal components
           and the cosine content.

           The cosine contents are stored as a ccs attribute of each covariance
           matrix object. The cosine content is a measure for random walk motion
           in the trajectory data. More background can be found in the work of
           Berk Hess:
           http://dx.doi.org/10.1103/PhysRevE.65.031910
           http://dx.doi.org/10.1103/PhysRevE.62.8438


           Argument:
             data  --  a new array with input data: rows correspond to time
                       frames and columns to coordinates/fields/...

           Optional argument
             output_prefix  --  A filename prefix for the tracks that will
                                contain the principal components. If not given,
                                the principal components are not written to
                                disk.

           Make sure the number of columns is always the same for one
           CovarianceBlocks instance.
        """
        def get_output_prefix():
            if output_prefix is None:
                return None
            else:
                if len(self.blocks) == 1:
                    return output_prefix
                else:
                    return "%s.%07i" % (output_prefix, self.proj_counter)

        if self.proj_counter == -1:
            self.proj_counter = 0
            self.blocks[0].init_proj(get_output_prefix())
        pos = 0
        while pos < len(data):
            #print "self.proj_counter", self.proj_counter
            #print "len(self.blocks)", len(self.blocks)
            #print "len(data)", len(data)
            #print "pos", pos
            #print "self.block_size", self.block_size
            #print
            if self.proj_counter >= len(self.blocks):
                return
                #raise ValueError("Can not project more data than was provided through add_data.")
            remaining = len(data) - pos
            todo = self.block_size - self.proj_done
            if remaining >= todo:
                self.blocks[self.proj_counter].data_proj(data[pos:pos+todo])
                #print self.blocks[self.proj_counter], "finish"
                self.blocks[self.proj_counter].finish_proj()
                pos += todo
                self.proj_done += todo
                self.proj_counter += 1
                if self.proj_counter < len(self.blocks):
                    self.proj_done = 0
                    self.blocks[self.proj_counter].init_proj(get_output_prefix())
            else: # todo > remaining
                self.blocks[self.proj_counter].data_proj(data[pos:])
                self.proj_done += remaining
                break

    def get_averages(self):
        """Compute the most relevant global properties over the blocks."""

        all_sigmas = numpy.array([block.sigmas for block in self.blocks])
        all_ccs = numpy.array([block.ccs for block in self.blocks])

        sigmas = all_sigmas.mean(axis=0)
        ccs = all_ccs.mean(axis=0)
        if len(self.blocks) > 1:
            sigmas_err = all_sigmas.std(axis=0)/len(self.blocks)
            ccs_err = all_ccs.std(axis=0)/len(self.blocks)
            overlap_multi = cov_overlap_multi([block.cov for block in self.blocks])
        else:
            sigmas_err = None
            ccs_err = None
            overlap_multi = None

        return sigmas, sigmas_err, ccs, ccs_err, overlap_multi


def cov_overlap(A, B):
    """Compute the overlap between two covariance matrices.

       A and B are two square matrices of the same size (numpy arrays). The
       return value is a scalar in the range [0,1]. When the result is zero,
       the matrices are each others opposites, note that this will never happen
       for covariance matrices because such matrices are positive definite. When
       the result is one, both matrices are identical.

       The exact formula is derived in the following paper:
       Hess, B. Physical Review E 2002, 65, 031910
       Cite this paper when you use this routine.
    """
    distance_max = numpy.sqrt(numpy.trace(A) + numpy.trace(B))
    evals_A, evecs_A = numpy.linalg.eigh(A)
    evals_B, evecs_B = numpy.linalg.eigh(B)
    tmp_A = evecs_A * numpy.sqrt(abs(evals_A))
    tmp_B = evecs_B * numpy.sqrt(abs(evals_B))
    root_A = numpy.dot(evecs_A.transpose(), tmp_A)
    root_B = numpy.dot(evecs_B.transpose(), tmp_B)
    distance = numpy.linalg.norm((root_A - root_B).ravel())
    return 1-distance/distance_max


def cov_overlap_multi(covs):
    """Compute the similarity for a list of covariance matrices.

       covs is a list of square matrices with at least two elements (numpy
       arrays). The similarity between the matrices is computed as an RMSD over
       RMSA ratio, where RSMD is the RMS with respect to the average of the
       square root of the matrices and RMSA is an RMS with respect to zero. (A
       stands for absolute.)
    """
    if len(covs) < 2:
        raise ValueError("At least two covariance matrices are expected")
    covs_sqrt = numpy.zeros((len(covs), covs[0].shape[0], covs[0].shape[1]), float)
    for i, cov in enumerate(covs):
        evals, evecs = numpy.linalg.eigh(cov)
        tmp = evecs * numpy.sqrt(abs(evals))
        covs_sqrt[i] = numpy.dot(evecs.transpose(), tmp)

    average = covs_sqrt.mean(axis=0)
    msd = ((covs_sqrt - average)**2).mean()
    msa = (covs_sqrt**2).mean()
    return numpy.sqrt(msd/msa)


pca_common_usage = """The following files are always written to the tracks database
${output_prefix}.cov        : a flattened covariance matrix (NxN elements)
${output_prefix}.evals      : the covariance eigen values
${output_prefix}.sigmas     : the square roots of the covariance eigen values
${output_prefix}.mode.${i}  : the covariance modes
${output_prefix}.ccs        : the cosine content of each mode
${output_prefix}.cosamp     : the amplitudes of the cosines of each mode

Optionally, the principal components (the time dependent amplitudes) can be
written to disk.
${output_prefix}.pc.${i}    : the principal components
"""


def pca_common_script(
    paths_in, dtype, sub, weights, correlation, reference, output_prefix,
    num_levels, dump_pcs, unit_name, unit
):
    """Shared code by tr-pca and tr-pca-geom.

       This function has only one purpose: bugs common to tr-pca and tr-pca-geom
       have to be fixed only once.
    """
    if num_levels < 2:
        raise ValueError("num_levels must be at least 2.")
    # call pca routine in tracks.api
    mtr = MultiTracksReader(paths_in, dtype, sub=sub)
    cm, overlap, cbs = pca_levels(mtr, num_levels, weights, correlation, reference)

    # dump some stuff to disk
    dump_track("%s.cov" % output_prefix, cm.cov.ravel())
    dump_track("%s.evals" % output_prefix, cm.evals)
    dump_track("%s.sigmas" % output_prefix, cm.sigmas)
    for i in xrange(len(cm.evals)):
        dump_track("%s.mode.%07i" % (output_prefix, i), cm.evecs[:,i])

    # compute the cosine contents and optionally write the principal components to disk
    mtr = MultiTracksReader(paths_in, dtype, sub=sub)
    for data in mtr.iter_buffers():
        data = data["data"]
        for level in xrange(num_levels):
            if level==0 and dump_pcs:
                cbs[level].project_data(data, output_prefix)
            else:
                cbs[level].project_data(data)
    dump_track("%s.ccs" % output_prefix, cm.ccs)
    dump_track("%s.cosamp" % output_prefix, cm.dot_data_cos/cm.sqnorms_cos)

    # Print some nice screen output with the most relevant results
    print "Overlap between the covariance of the first and the second half of the trajectory:"
    print "  %.2f %%" % (overlap*100)
    print
    white = (" "*len(unit_name))
    for level in xrange(num_levels):
        sigmas, sigmas_err, ccs, ccs_err, overlap_multi = cbs[level].get_averages()
        print "Level %i: averages over blocks with size total/%i=%i" % (
            level, 2**level, mtr.shortest/(2**level)
        )
        if level > 0:
            print "RMSD/RMSA of the covariance matrices of each block: %.2f %%" % (overlap_multi*100)
            print
        if level == 0:
            print "            Sigma [%s]    Cosine content  " % unit_name
            print "           ------------%s  ---------------- " % white
            for i in xrange(len(cm.sigmas)):
                print " %7i    %9.3e%s    %8.2f %%" % (
                    i, sigmas[i]/unit, white, ccs[i]*100
                )
        else:
            print "               Sigma [%s]                   Cosine content  " % unit_name
            print "          ----------------------%s      ---------------------- " % white
            for i in xrange(len(cm.sigmas)):
                print " %7i   %9.3e +- %9.3e%s    %8.2f %% +- %6.2f %%" % (
                    i, sigmas[i]/unit, sigmas_err[i]/unit, white,
                    ccs[i]*100, ccs_err[i]*100
                )
        print
        print
        if level > 0:
            dump_track("%s.ccs.%07i" % (output_prefix, level), ccs)
            dump_track("%s.ccs_err.%07i" % (output_prefix, level), ccs_err)
            dump_track("%s.sigmas.%07i" % (output_prefix, level), sigmas)
            dump_track("%s.sigmas_err.%07i" % (output_prefix, level), sigmas_err)

    return cm.mean

