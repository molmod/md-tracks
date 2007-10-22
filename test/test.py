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

import os, sys

orig_dir = os.getcwd()
scripts_dir = os.path.join(os.path.dirname(os.getcwd()), "scripts")
lib_dir = os.path.join(os.path.dirname(os.getcwd()), "src")
tmp_dir = os.path.join(os.getcwd(), "tmp")
input_dir = os.path.join(os.getcwd(), "input")
output_dir = os.path.join(os.getcwd(), "output")
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

sys.path.insert(0, lib_dir)

import unittest, os, shutil, glob
from ccio.psf import PSFFile
from ccio.xyz import XYZReader, XYZFile
from molmod.units import angstrom, fs
from molmod.constants import lightspeed
from molmod.data import periodic
from tracks import load_track, dist_track, bend_track, dihed_track, dump_track, PSFFilter
import numpy


class CommandsTestCase(unittest.TestCase):
    def setUp(self):
        os.makedirs(tmp_dir)
        os.chdir(tmp_dir)

    def tearDown(self):
        os.chdir(orig_dir)
        shutil.rmtree(tmp_dir)

    def from_xyz(self, case, middle_word, *extra_args):
        self.execute("tr-from-xyz", [
            os.path.join(input_dir, case, "md-%s-1.xyz" % middle_word),
            middle_word
        ] + list(extra_args))

    def from_cp2k_ener(self, case, *extra_args):
        self.execute("tr-from-cp2k-ener", [
            os.path.join(input_dir, case, "md-1.ener"),
        ] + list(extra_args))

    def from_cpmd_traj(self, filename="cpmd_traj", *extra_args):
        self.execute("tr-from-cpmd-traj", [
            os.path.join(input_dir, filename),
        ] + list(extra_args))

    def execute(self, command, args, verbose=False, stdin=None):
        from subprocess import Popen, PIPE, STDOUT
        p = Popen(
            ["/usr/bin/python", os.path.join(scripts_dir, command)] + args,
            stdin=PIPE, stdout=PIPE, stderr=STDOUT, env={"PYTHONPATH": lib_dir},
        )
        if stdin is not None:
            for line in stdin:
                print >> p.stdin, line
        p.stdin.close()
        output = list(p.stdout)
        retcode = p.wait()
        self.assert_(retcode==0, "Something went wrong with this command:\n%s %s\n. The output is:\n%s" % (command, " ".join(args), "".join(output)))
        if verbose:
            print "".join(output)
        return output

    def test_from_xyz(self):
        # Load the xyz file
        self.from_xyz("thf01", "pos")
        # Test some values
        tmp = load_track("tracks/atom.pos.0000000.x")
        self.assertAlmostEqual(tmp[0]/angstrom, 1.160855, 5)
        self.assertAlmostEqual(tmp[1]/angstrom, 1.1814236022, 5)
        self.assertAlmostEqual(tmp[-1]/angstrom, -0.9141294178, 5)
        tmp = load_track("tracks/atom.pos.0000012.z")
        self.assertAlmostEqual(tmp[0]/angstrom, 1.0237910000, 5)
        self.assertAlmostEqual(tmp[1]/angstrom, 1.0193867160, 5)
        self.assertAlmostEqual(tmp[-1]/angstrom, -0.5994763399, 5)
        # Load
        self.from_xyz("thf01", "pos", "-s20:601:5")
        # Test some values
        tmp = load_track("tracks/atom.pos.0000000.x")
        self.assertAlmostEqual(tmp[0]/angstrom, 1.1643775386, 5)
        self.assertAlmostEqual(tmp[1]/angstrom, 1.1186662255, 5)
        self.assertAlmostEqual(tmp[-1]/angstrom, 0.3461355118, 5)
        tmp = load_track("tracks/atom.pos.0000012.z")
        self.assertAlmostEqual(tmp[0]/angstrom, 1.4974560873, 5)
        self.assertAlmostEqual(tmp[1]/angstrom, 1.7383838088, 5)
        self.assertAlmostEqual(tmp[-1]/angstrom, -0.6220795393, 5)
        # Load the xyz file
        self.from_xyz("thf01", "vel", "-u1")
        # Test some values
        tmp = load_track("tracks/atom.vel.0000000.x")
        self.assertAlmostEqual(tmp[0], 0.0002092059, 5)
        self.assertAlmostEqual(tmp[1], 0.0001447176, 5)
        self.assertAlmostEqual(tmp[-1], -0.0001092007, 5)
        tmp = load_track("tracks/atom.vel.0000012.z")
        self.assertAlmostEqual(tmp[0], -0.0003909859, 5)
        self.assertAlmostEqual(tmp[1], 0.0004487963, 5)
        self.assertAlmostEqual(tmp[-1], 0.0000859001, 5)
        # Load the xyz file
        self.from_xyz("thf01", "vel", "-u1", "-s20:601:5")
        # Test some values
        tmp = load_track("tracks/atom.vel.0000000.x")
        self.assertAlmostEqual(tmp[0], 0.0000503137, 5)
        self.assertAlmostEqual(tmp[1], -0.0000955072, 5)
        self.assertAlmostEqual(tmp[-1], -0.0000947954, 5)
        tmp = load_track("tracks/atom.vel.0000012.z")
        self.assertAlmostEqual(tmp[0], 0.0001216775, 5)
        self.assertAlmostEqual(tmp[1], -0.0001251857, 5)
        self.assertAlmostEqual(tmp[-1], 0.0007943491, 5)

    def test_from_cp2k_ener(self):
        # Load the energy file
        self.from_cp2k_ener("thf01")
        # Test some values
        tmp = load_track("tracks/step")
        self.assertEqual(tmp[0], 0)
        self.assertEqual(tmp[1], 5)
        self.assertEqual(tmp[-1], 5000)
        tmp = load_track("tracks/time")
        self.assertAlmostEqual(tmp[0]/fs, 0.0, 5)
        self.assertAlmostEqual(tmp[1]/fs, 5.0, 5)
        self.assertAlmostEqual(tmp[-1]/fs, 5000.0, 5)
        tmp = load_track("tracks/kinetic_energy")
        self.assertAlmostEqual(tmp[0], 0.015675735, 5)
        self.assertAlmostEqual(tmp[1], 0.008711175, 5)
        self.assertAlmostEqual(tmp[-1], 0.011256845, 5)
        tmp = load_track("tracks/temperature")
        self.assertAlmostEqual(tmp[0], 300.000000000, 5)
        self.assertAlmostEqual(tmp[1], 166.713237534, 5)
        self.assertAlmostEqual(tmp[-1], 215.431909415, 5)
        tmp = load_track("tracks/potential_energy")
        self.assertAlmostEqual(tmp[0], 0.029894724, 5)
        self.assertAlmostEqual(tmp[1], 0.036975396, 5)
        self.assertAlmostEqual(tmp[-1], 0.034406422, 5)
        tmp = load_track("tracks/total_energy")
        self.assertAlmostEqual(tmp[0], 0.045570459, 5)
        self.assertAlmostEqual(tmp[1], 0.045686571, 5)
        self.assertAlmostEqual(tmp[-1], 0.045663267, 5)
        # Load the energy file
        self.from_cp2k_ener("thf01", "-s20:601:5")
        # Test some values
        tmp = load_track("tracks/step")
        self.assertEqual(tmp[0], 100)
        self.assertEqual(tmp[1], 125)
        self.assertEqual(tmp[-1], 3000)
        tmp = load_track("tracks/time")
        self.assertAlmostEqual(tmp[0]/fs, 100.0, 5)
        self.assertAlmostEqual(tmp[1]/fs, 125.0, 5)
        self.assertAlmostEqual(tmp[-1]/fs, 3000.0, 5)
        tmp = load_track("tracks/kinetic_energy")
        self.assertAlmostEqual(tmp[0], 0.007844166, 5)
        self.assertAlmostEqual(tmp[1], 0.008257514, 5)
        self.assertAlmostEqual(tmp[-1], 0.011402720, 5)
        tmp = load_track("tracks/temperature")
        self.assertAlmostEqual(tmp[0], 150.120545659, 5)
        self.assertAlmostEqual(tmp[1], 158.031128585, 5)
        self.assertAlmostEqual(tmp[-1], 218.223631032, 5)
        tmp = load_track("tracks/potential_energy")
        self.assertAlmostEqual(tmp[0], 0.037852247, 5)
        self.assertAlmostEqual(tmp[1], 0.037445922, 5)
        self.assertAlmostEqual(tmp[-1], 0.034186532, 5)
        tmp = load_track("tracks/total_energy")
        self.assertAlmostEqual(tmp[0], 0.045696413, 5)
        self.assertAlmostEqual(tmp[1], 0.045703436, 5)
        self.assertAlmostEqual(tmp[-1], 0.045589252, 5)

    def test_to_xyz(self):
        self.from_xyz("thf01", "pos")
        self.execute("tr-to-xyz", [os.path.join(input_dir, "thf01/init.xyz"), "tracks/atom.pos", "test.pos.xyz"])
        xyz_reader_orig = XYZReader(os.path.join(input_dir, "thf01/md-pos-1.xyz"))
        xyz_reader_copy = XYZReader("test.pos.xyz")
        self.assert_(xyz_reader_orig.symbols == xyz_reader_copy.symbols)
        for (title_orig, coordinates_orig), (tile_copy, coordinates_copy) in zip(xyz_reader_orig, xyz_reader_copy):
            self.assert_(abs(coordinates_orig - coordinates_copy).max() < 1e-7)

    def test_read_write_slice_length(self):
        self.from_cp2k_ener("water32")
        # tr-read
        lines = self.execute("tr-read", ["tracks/temperature"])
        tmp1 = numpy.array([float(line) for line in lines], float)
        tmp2 = load_track("tracks/temperature")
        self.assert_((tmp1==tmp2).all())
        self.assert_(len(lines)==201)
        # tr-length
        length = int("".join(self.execute("tr-length", ["tracks/temperature"])))
        self.assert_(length==201)
        # tr-slice
        self.execute("tr-slice", ["tracks/temperature", "20:80:5", "tracks/temperature_sliced"])
        length = int("".join(self.execute("tr-length", ["tracks/temperature_sliced"])))
        self.assert_(length==12)
        t = load_track("tracks/temperature")
        ts = load_track("tracks/temperature_sliced")
        self.assert_((t[20:80:5] == ts).all())
        # tr-write
        tmp1 = numpy.arange(101, dtype=float)
        lines = [str(val) for val in tmp1]
        self.execute("tr-write", ["tracks/tmp"], stdin=lines)
        tmp2 = load_track("tracks/tmp")
        self.assert_((tmp1==tmp2).all())

    def test_ac(self):
        self.from_xyz("thf01", "vel", "-u1")
        self.from_cp2k_ener("thf01")
        self.execute("tr-ac", glob.glob("tracks/atom.vel.*") + ["--tau=200*fs", "5.0*fs", "tracks/vac_a1"])
        self.execute("tr-ac", glob.glob("tracks/atom.vel.*") + ["--tau=200*fs", "tracks/time", "tracks/vac_a2"])
        length = int("".join(self.execute("tr-length", ["tracks/vac_a2"])))
        self.assert_(length==501)
        self.execute("tr-slice", ["tracks/time", ":%i:" % length, "tracks/time_sliced"])
        self.execute("tr-plot", [
            "tracks/time_sliced,tracks/vac_a1,tracks/vac_a1.error",
            "-x", "delta t", "-y", "VAC", "-t", "Velocity autocorrelation function (thf01)",
            "--xunit=ps", os.path.join(output_dir, "ac_vac_a1")
        ])
        self.execute("tr-plot", [
            "tracks/time_sliced,tracks/vac_a1.normalized,tracks/vac_a1.normalized.error",
            "--ylim=-1,1", "-x", "delta t", "-y", "VAC", "-t", "Normalized velocity autocorrelation function (thf01)",
            "--xunit=ps", os.path.join(output_dir, "ac_vac_a1.normalized")
        ])
        tmp1 = load_track("tracks/vac_a1")
        tmp2 = load_track("tracks/vac_a2")
        self.assert_((tmp1==tmp2).all())
        tmp1 = load_track("tracks/vac_a1.normalized")
        tmp2 = load_track("tracks/vac_a2.normalized")
        self.assert_((tmp1==tmp2).all())
        self.assert_(tmp1[0]==1.0)
        tmp1 = load_track("tracks/vac_a1.error")
        tmp2 = load_track("tracks/vac_a2.error")
        self.assert_((tmp1==tmp2).all())
        self.assert_((tmp1/tmp1[0] == 1.0).all())
        tmp1 = load_track("tracks/vac_a1.normalized.error")
        tmp2 = load_track("tracks/vac_a2.normalized.error")
        self.assert_((tmp1==tmp2).all())
        self.assert_((tmp1/tmp1[0] == 1.0).all())

        self.execute("tr-ac", glob.glob("tracks/atom.vel.*") + ["--tau=200*fs", "-m 4000*fs", "5.0*fs", "tracks/vac_b1"])
        self.execute("tr-ac", glob.glob("tracks/atom.vel.*") + ["--tau=200*fs", "-m 4000*fs", "tracks/time", "tracks/vac_b2"])
        length = int("".join(self.execute("tr-length", ["tracks/vac_b2"])))
        self.assert_(length==801)
        self.execute("tr-slice", ["tracks/time", ":%i:" % length, "tracks/time_sliced"])
        self.execute("tr-plot", [
            "tracks/time_sliced,tracks/vac_b1,tracks/vac_b1.error",
            "-x", "delta t", "-y", "VAC", "-t", "Velocity autocorrelation function (thf01)",
            "--xunit=ps", os.path.join(output_dir, "ac_vac_b1")
        ])
        self.execute("tr-plot", [
            "tracks/time_sliced,tracks/vac_b1.normalized,tracks/vac_b1.normalized.error",
            "--ylim=-1,1", "-x", "delta t", "-y", "VAC", "-t", "Normalized velocity autocorrelation function (thf01)",
            "--xunit=ps", os.path.join(output_dir, "ac_vac_b1.normalized")
        ])
        tmp1 = load_track("tracks/vac_b1")
        tmp2 = load_track("tracks/vac_b2")
        self.assert_((tmp1==tmp2).all())
        tmp1 = load_track("tracks/vac_b1.normalized")
        tmp2 = load_track("tracks/vac_b2.normalized")
        self.assert_((tmp1==tmp2).all())
        self.assert_(tmp1[0]==1.0)
        tmp1 = load_track("tracks/vac_b1.error")
        tmp2 = load_track("tracks/vac_b2.error")
        self.assert_((tmp1==tmp2).all())
        self.assert_((tmp1/tmp1[0] == 1.0).all())
        tmp1 = load_track("tracks/vac_b1.normalized.error")
        tmp2 = load_track("tracks/vac_b2.normalized.error")
        self.assert_((tmp1==tmp2).all())
        self.assert_((tmp1/tmp1[0] == 1.0).all())

    def test_integrate(self):
        self.from_xyz("thf01", "vel", "-u1")
        self.from_cp2k_ener("thf01")
        self.execute("tr-ac", glob.glob("tracks/atom.vel.*") + ["--tau=200*fs", "-m3000*fs", "tracks/time", "tracks/vac"])
        self.execute("tr-integrate", ["tracks/vac.normalized", "tracks/time"])
        length = int("".join(self.execute("tr-length", ["tracks/vac.normalized"])))
        self.assert_(length==601)
        length = int("".join(self.execute("tr-length", ["tracks/vac.normalized.error"])))
        self.assert_(length==601)
        length = int("".join(self.execute("tr-length", ["tracks/vac.normalized.int"])))
        self.assert_(length==601)
        length = int("".join(self.execute("tr-length", ["tracks/vac.normalized.int.error"])))
        self.assert_(length==601)
        self.execute("tr-slice", ["tracks/time", ":%i:" % length, "tracks/time_sliced"])
        self.execute("tr-plot", [
            "tracks/time_sliced,tracks/vac.normalized,tracks/vac.normalized.error",
            "--ylim=-1,1", "-x", "delta t", "-y", "VAC", "-t", "Normalized velocity autocorrelation function (thf01)",
            "--xunit=ps", os.path.join(output_dir, "integrate_vac.normalized")
        ])
        self.execute("tr-plot", [
            "tracks/time_sliced,tracks/vac.normalized.int,tracks/vac.normalized.int.error",
            "-x", "delta t", "-y", "Int(VAC)", "-t", "Integral of the normalized velocity autocorrelation function (thf01)",
            "--xunit=ps", "--yunit=fs", os.path.join(output_dir, "integrate_vac.normalized.int")
        ])

    def test_rfft_irfft(self):
        self.from_xyz("thf01", "vel", "-s:1000:", "-u1") # irrft always results in an even number of datapoints
        self.from_cp2k_ener("thf01")
        self.execute("tr-rfft", glob.glob("tracks/atom.vel.*.?"))
        self.execute("tr-irfft", glob.glob("tracks/atom.vel.*.rfft"))
        self.assert_(len(glob.glob("tracks/atom.vel.*.?")) == len(glob.glob("tracks/*.rfft")))
        self.assert_(len(glob.glob("tracks/atom.vel.*.?")) == len(glob.glob("tracks/*.rfft.irfft")))
        for filename in glob.glob("tracks/atom.vel.*.?"):
            other_filename = "%s.rfft.irfft" % filename
            tmp1 = load_track(filename)
            tmp2 = load_track(other_filename)
            self.assert_(tmp1.shape == tmp2.shape)
            self.assert_(abs(tmp1-tmp2).max() < 1e-7)

    def test_make_spectrum(self):
        self.from_xyz("thf01", "vel", "-u1")
        self.from_cp2k_ener("thf01")
        self.execute("tr-rfft", glob.glob("tracks/atom.vel.*") + ["-whamming", "-d 5"])
        self.execute("tr-make-spectrum", glob.glob("tracks/atom.vel.*.rfft") + ["tracks/spectrum"])
        self.assert_(len(load_track("tracks/spectrum"))==101)
        self.execute("tr-wavenumber-axis", ["tracks/spectrum", "1000*fs", "tracks/wavenumbers"])
        self.execute("tr-freq-axis", ["tracks/spectrum", "1000*fs", "tracks/freqs"])
        wavenumbers = load_track("tracks/wavenumbers")
        freqs = load_track("tracks/freqs")
        self.assert_(wavenumbers[0] == 0.0)
        self.assert_(freqs[0] == 0.0)
        self.assert_((abs(freqs[1:]/wavenumbers[1:]-lightspeed)/lightspeed).max() < 1e-7)
        self.execute("tr-plot", [
            "--xlabel=Wavenumber", "-s1::", "--ylabel=Amplitude", "--xunit=1/cm",
            "tracks/wavenumbers,tracks/spectrum",
            os.path.join(output_dir, "make_spectrum_wavenumbers")]
        )
        self.execute("tr-plot", [
            "--xlabel=Frequency", "-s1::", "--ylabel=Amplitude", "--xunit=1/fs",
            "tracks/freqs,tracks/spectrum",
            os.path.join(output_dir, "make_spectrum_freqs")]
        )
        self.execute("tr-plot", [
            "--xlabel=Time", "-s1::", "--ylabel=Amplitude", "--xunit=fs", "--xinv",
            "tracks/freqs,tracks/spectrum",
            os.path.join(output_dir, "make_spectrum_freqs_inv")]
        )

    def test_fit_peaks(self):
        self.from_xyz("thf01", "vel", "-u1")
        self.from_cp2k_ener("thf01")
        self.execute("tr-rfft", glob.glob("tracks/atom.vel.*") + ["-whamming", "-d 5"])
        self.execute("tr-make-spectrum", glob.glob("tracks/atom.vel.*.rfft") + ["tracks/spectrum"])
        self.execute("tr-wavenumber-axis", ["tracks/spectrum", "1000*fs", "tracks/wavenumbers"])
        output = self.execute("tr-fit-peaks", [
            "tracks/wavenumbers", "tracks/spectrum", "1750", "2250",
            "0.0001:2050.0:50.0:0.01", "--dump-model=tracks/model", #"--no-fit",
        ])
        f = file(os.path.join(output_dir, "fit_peaks.out"), "w")
        f.writelines(output)
        f.close()
        self.execute("tr-plot", [
            "--xlabel=Wavenumber", "-s3::", "--ylabel=Amplitude",
            "--xunit=1/cm", "tracks/wavenumbers,tracks/spectrum",
            "tracks/wavenumbers,tracks/model",
            os.path.join(output_dir, "fit_peaks_spectrum"),
        ])

    def test_freq_axis(self):
        self.from_cp2k_ener("thf01")
        time = load_track("tracks/time")
        self.execute("tr-rfft", ["tracks/temperature"])
        self.execute("tr-make-spectrum", ["tracks/temperature.rfft", "tracks/spectrum"])
        spectrum = load_track("tracks/spectrum")
        self.execute("tr-freq-axis", ["tracks/spectrum", "tracks/time", "tracks/freqs"])
        freqs = load_track("tracks/freqs")
        self.assert_(len(freqs)==len(spectrum))
        self.assert_(freqs[0]==0)
        self.assert_(freqs[1]==1/time[-1])
        self.assert_(freqs[-1]==1/(2*(time[1]-time[0])))
        self.execute("tr-freq-axis", ["tracks/spectrum", "5000*fs", "tracks/freqs"])
        freqs_check = load_track("tracks/freqs")
        self.assert_((freqs==freqs_check).all())
        self.execute("tr-freq-axis", ["501", "tracks/time", "tracks/freqs"])
        freqs_check = load_track("tracks/freqs")
        self.assert_((freqs==freqs_check).all())
        self.execute("tr-freq-axis", ["501", "5000*fs", "tracks/freqs"])
        freqs_check = load_track("tracks/freqs")
        self.assert_((freqs==freqs_check).all())

    def test_wavenumber_axis(self):
        self.from_cp2k_ener("thf01")
        time = load_track("tracks/time")
        self.execute("tr-rfft", ["tracks/temperature"])
        self.execute("tr-make-spectrum", ["tracks/temperature.rfft", "tracks/spectrum"])
        spectrum = load_track("tracks/spectrum")
        self.execute("tr-wavenumber-axis", ["tracks/spectrum", "tracks/time", "tracks/wavenumbers"])
        wavenumbers = load_track("tracks/wavenumbers")
        self.assert_(len(wavenumbers)==len(spectrum))
        self.assert_(wavenumbers[0]==0)
        self.assert_(wavenumbers[1]==1/time[-1]/lightspeed)
        self.assert_(wavenumbers[-1]==1/(2*(time[1]-time[0]))/lightspeed)
        self.execute("tr-wavenumber-axis", ["tracks/spectrum", "5000*fs", "tracks/wavenumbers"])
        wavenumbers_check = load_track("tracks/wavenumbers")
        self.assert_((wavenumbers==wavenumbers_check).all())
        self.execute("tr-wavenumber-axis", ["501", "tracks/time", "tracks/wavenumbers"])
        wavenumbers_check = load_track("tracks/wavenumbers")
        self.assert_((wavenumbers==wavenumbers_check).all())
        self.execute("tr-wavenumber-axis", ["501", "5000*fs", "tracks/wavenumbers"])
        wavenumbers_check = load_track("tracks/wavenumbers")
        self.assert_((wavenumbers==wavenumbers_check).all())

    def test_ic_dist(self):
        self.from_xyz("thf01", "pos")
        self.execute("tr-ic-dist", ["tracks/atom.pos", "1", "2"])
        dists = load_track("tracks/atom.pos.dist.0000001.0000002")
        self.assertAlmostEqual(dists[0], 4.25631, 4)
        self.assertAlmostEqual(dists[1], 4.28458, 4)
        self.assertAlmostEqual(dists[-1], 4.26709, 4)
        self.execute("tr-ic-dist", ["-s20:601:5", "tracks/atom.pos", "1", "2"])
        dists = load_track("tracks/atom.pos.dist.0000001.0000002")
        self.assertAlmostEqual(dists[0], 4.32567, 4)
        self.assertAlmostEqual(dists[1], 4.41805, 4)
        self.assertAlmostEqual(dists[-1], 4.35014, 4)

    def test_ic_bend(self):
        self.from_xyz("thf01", "pos")
        self.execute("tr-ic-bend", ["tracks/atom.pos", "1", "0", "2"])
        bends = load_track("tracks/atom.pos.bend.0000001.0000000.0000002")
        self.assertAlmostEqual(bends[0]*180/numpy.pi, 105.426, 3)
        self.assertAlmostEqual(bends[1]*180/numpy.pi, 102.286, 3)
        self.assertAlmostEqual(bends[-1]*180/numpy.pi, 107.284, 3)
        self.execute("tr-ic-bend", ["-s20:601:5", "tracks/atom.pos", "1", "0", "2"])
        bends = load_track("tracks/atom.pos.bend.0000001.0000000.0000002")
        self.assertAlmostEqual(bends[0]*180/numpy.pi, 104.739, 3)
        self.assertAlmostEqual(bends[1]*180/numpy.pi, 108.972, 3)
        self.assertAlmostEqual(bends[-1]*180/numpy.pi, 107.277, 3)

    def test_ic_dihed(self):
        self.from_xyz("thf01", "pos")
        self.execute("tr-ic-dihed", ["tracks/atom.pos", "2", "3", "4", "1"])
        bends = load_track("tracks/atom.pos.dihed.0000002.0000003.0000004.0000001")
        self.assertAlmostEqual(bends[0]*180/numpy.pi, 0.000, 3)
        self.assertAlmostEqual(bends[1]*180/numpy.pi, -1.919, 3)
        self.assertAlmostEqual(bends[-1]*180/numpy.pi, 21.320, 3)
        self.execute("tr-ic-dihed", ["-s20:601:5", "tracks/atom.pos", "2", "3", "4", "1"])
        bends = load_track("tracks/atom.pos.dihed.0000002.0000003.0000004.0000001")
        self.assertAlmostEqual(bends[0]*180/numpy.pi, -15.209, 3)
        self.assertAlmostEqual(bends[1]*180/numpy.pi, -16.602, 3)
        self.assertAlmostEqual(bends[-1]*180/numpy.pi, -31.306, 3)

    def test_ic_psf(self):
        self.from_xyz("thf01", "pos")
        self.execute("tr-ic-psf", ["tracks/atom.pos", os.path.join(input_dir, "thf01/init.psf")])
        # bond
        bond_filenames = glob.glob("tracks/atom.pos.dist*")
        self.assert_(len(bond_filenames)==13)
        for bond_filename in bond_filenames:
            bond = load_track(bond_filename)
            index1, index2 = [int(word) for word in bond_filename.split(".")[-2:]]
            bond_check = dist_track(index1, index2, "tracks/atom.pos", slice(None))
            self.assert_((bond==bond_check).all())
        # bend
        bend_filenames = glob.glob("tracks/atom.pos.bend*")
        self.assert_(len(bend_filenames)==25)
        for bend_filename in bend_filenames:
            bend = load_track(bend_filename)
            index1, index2, index3 = [int(word) for word in bend_filename.split(".")[-3:]]
            bend_check = bend_track(index1, index2, index3, "tracks/atom.pos", slice(None))
            self.assert_((bend==bend_check).all())
        # dihed
        dihed_filenames = glob.glob("tracks/atom.pos.dihed*")
        self.assert_(len(dihed_filenames)==33)
        for dihed_filename in dihed_filenames:
            dihed = load_track(dihed_filename)
            index1, index2, index3, index4 = [int(word) for word in dihed_filename.split(".")[-4:]]
            dihed_check = dihed_track(index1, index2, index3, index4, "tracks/atom.pos", slice(None))
            self.assert_((dihed==dihed_check).all())

    def test_mean_std(self):
        self.from_xyz("thf01", "vel", "-u1")
        self.from_cp2k_ener("thf01")
        molecule = XYZFile(os.path.join(input_dir, "thf01/init.xyz")).get_molecule()
        # A1) compute the kinetic energy per atom:
        for index, number in enumerate(molecule.numbers):
            vx = load_track("tracks/atom.vel.%07i.x" % index)
            vy = load_track("tracks/atom.vel.%07i.y" % index)
            vz = load_track("tracks/atom.vel.%07i.z" % index)
            ekin = 0.5*periodic[number].mass*(vx*vx+vy*vy+vz*vz)
            dump_track("tracks/atom.ekin.%07i" % index, ekin)
        # A2) average and compare to the kinetic energy from the energy file
        self.execute("tr-mean-std", glob.glob("tracks/atom.ekin.*") + ["tracks/atom.ekin"])
        ekin_mean = load_track("tracks/atom.ekin.mean")
        ekin = load_track("tracks/kinetic_energy")
        self.assert_(abs(ekin/13-ekin_mean).max()/ekin_mean.max() < 1e-6)
        # B) Verify the relation between std and error
        self.execute("tr-mean-std", glob.glob("tracks/atom.vel.*") + ["tracks/atom.vel"])
        vel_error = load_track("tracks/atom.vel.error")
        vel_std = load_track("tracks/atom.vel.std")
        self.assert_(abs(vel_std-vel_error*numpy.sqrt(3*13)).max()/abs(vel_std).max() < 1e-5)

    def test_blav(self):
        self.from_cp2k_ener("thf01")
        self.execute("tr-blav", ["tracks/temperature", "tracks/time", "-b10", "-tfs"])
        self.execute("tr-blav", [
            "tracks/temperature", "tracks/time", "-b5", "-tfs",
            "--plot_error=%s" % os.path.join(output_dir, "blav_error"),
            "--plot_ctime=%s" % os.path.join(output_dir, "blav_ctime"),
        ])

    def test_split_com(self):
        self.from_xyz("water32", "vel", "-u1")
        self.from_cp2k_ener("water32")
        self.execute("tr-split-com", ["tracks/atom.vel", "vel", os.path.join(input_dir, "water32/init.psf")])
        psf = PSFFile(os.path.join(input_dir, "water32/init.psf"))
        # check that the coms have in total no translational kinetic energy
        for c in 'xyz':
            self.execute("tr-mean-std", glob.glob("tracks/com.vel.*.%s" % c) + ["tracks/com.vel.%s" % c])
            vmean = load_track("tracks/com.vel.%s.mean" % c)
            self.assert_(abs(vmean).max() < 1e-10)
        # check that the weighted sum of the relative velocities is always zero, per molecule
        for m_index in xrange(32):
            for c in 'xyz':
                tmp = 0
                for a_index in (psf.molecules==m_index).nonzero()[0]:
                    mass = periodic[psf.numbers[a_index]].mass
                    tmp += mass*load_track("tracks/rel.vel.%07i.%s" % (a_index, c))
                self.assert_(abs(tmp).max() < 1e-14)
        # compute the total kinetic energy from the com and the rel contributions
        # and compare it to tracks/kinetic_energy
        ekin = 0
        mass_water = periodic[1].mass*2+periodic[8].mass
        for m_index in xrange(32):
            for c in 'xyz':
                v = load_track("tracks/com.vel.%07i.%s" % (m_index, c))
                ekin += 0.5*mass_water*v*v
                for a_index in (psf.molecules==m_index).nonzero()[0]:
                    mass = periodic[psf.numbers[a_index]].mass
                    v = load_track("tracks/rel.vel.%07i.%s" % (a_index, c))
                    ekin += 0.5*mass*v*v
        ekin_check = load_track("tracks/kinetic_energy")
        self.assert_(abs(ekin/ekin_check - 1).max() < 1e-6)

    def test_filter(self):
        verbose = False
        if verbose: print
        def check_filter(case, kind, expression, expected):
            arguments = [os.path.join(input_dir, "%s/init.psf" % case), kind, expression]
            result = self.execute("tr-filter", arguments)[0].strip()
            self.assertEqual(result, expected)
            if verbose: print "%s  |  %s  |  %s   =>   %s" % (case, kind, expression, result)
            arguments.append("--prefix=test")
            result = self.execute("tr-filter", arguments)[0].strip()
            if verbose: print "%s  |  %s  |  %s   =>   %s" % (case, kind, expression, result)
            arguments.append("--xyz")
            result = self.execute("tr-filter", arguments)[0].strip()
            if verbose: print "%s  |  %s  |  %s   =>   %s" % (case, kind, expression, result)

        check_filter('thf01', 'at', 'a.symbol=="c"', '1,2,3,4')
        check_filter('thf01', 'at', 'a.nsymbols=="c,c,h,h"', '3,4')
        check_filter('thf01', 'at', 'a.nsymbols=="o,c,h,h"', '1,2')
        check_filter('thf01', 'at', 'a.nsymbols=="c_2,h_2"', '3,4')
        check_filter('thf01', 'at', 'a.nnumbers=="6_2,1_2"', '3,4')
        check_filter('thf01', 'at', 'a.nnumbers=="6,8,1,1"', '1,2')
        check_filter('thf01', 'at', 'a.nnumbers=="6,6"', '0')
        check_filter('thf01', 'at', 'a.number=="6"', '')
        check_filter('thf01', 'at', 'a.number==6', '1,2,3,4')
        check_filter('thf01', 'at', 'a.index==6', '6')
        check_filter('thf01', 'at', '0 in a.nindexes', '1,2')
        check_filter('thf01', 'at', 'a.m.index==0', '0,1,2,3,4,5,6,7,8,9,10,11,12')
        check_filter('thf01', 'mol', 'a.index==6', '0')
        check_filter('thf01', 'mol', 'a.index=="6"', '')
        check_filter('thf01', 'mol', 'm.index==0', '0')
        check_filter('thf01', 'mol', 'm.composition=="C_4,O,H_8"', '0')
        check_filter('thf01', 'mol', 'm.composition=="C_3,O,H_8"', '')
        check_filter('water32', 'mol', 'm.composition=="H_2,O"', ",".join(str(val) for val in xrange(32)))
        check_filter('water32', 'mol', 'm.index==5', '5')
        check_filter('water32', 'mol', 'a.index==6', '2')



class PSFFilterTestCase(unittest.TestCase):
    def test_filter(self):
        psf = PSFFile(os.path.join(input_dir, "water32/init.psf"))

        psf_filter = PSFFilter(psf,[1,2,3,54],[])
        self.assert_(psf_filter(4, 3))
        self.assert_(not psf_filter(4, 7, 13))

        psf_filter = PSFFilter(psf,[],[1])
        self.assert_(psf_filter(3, 55, 77))
        self.assert_(not psf_filter(88, 22, 0))

        psf_filter = PSFFilter(psf,[1,2,10,54],[1])
        self.assert_(psf_filter(2, 88, 10))
        self.assert_(psf_filter(10, 0, 54))
        self.assert_(not psf_filter(30, 31, 0))


unittest.main()


