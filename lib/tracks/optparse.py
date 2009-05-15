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



def add_quiet_option(parser):
    parser.add_option(
        "-q", "--quiet", action="store_false", dest="verbose", default=True,
        help="Don't print any output."
    )

def add_slice_option(parser):
    parser.add_option(
        "-s", "--slice", default="::",
        help="Subsample the (time dependent) input tracks with the given slice "
             "start:stop:step where start, stop and step must be integers or "
             "can be omitted. The slice interpretation is pythonic. "
             "[default=%default]",
    )

def add_append_option(parser):
    parser.add_option(
        "--append", action="store_false", dest="clear", default=True,
        help="Append to existing tracks if possible."
    )

def add_cell_option(parser):
    parser.add_option(
        "-c", "--cell", dest="unit_cell_str", default=None,
        help="Take into account the given periodic boundary conditions. "
             "The unit cell parameters. Several formats are supported. "
             "(i) 'a,' A cubic unit cell with ridge a. "
             "(ii) 'a,b,c' The parameters of an orthorhombic cell. "
             "(iii) 'a,b,c,alpha,beta,gamma' The parameters for a triclinic cell. "
             "(iv) 'ax,ay,az,bx,by,bz,cx,cy,cz' The cartesian parameters for a triclinic cell. "
             "(v) 'cell_prefix' A track prefix can be used for a time dependent unit cell. "
             "The presence of comma's is used to differentiate between all "
             "the possibilities."
    )

def add_cos_option(parser):
    parser.add_option(
        "--cos", action="store_true", default=False,
        help="Compute the cosine instead of the angle."
    )

def add_filter_atoms_option(parser):
    parser.add_option(
        "-a", "--filter-atoms",
        help="Only consider the atoms listed in FILTER_ATOMS. FILTER_ATOMS is a "
             "comma separated list of of integers. Counting starts at zero.",
    )

def add_ic_project_option(parser, name):
    parser.add_option(
        "-p", "--project", action="store_true", default=False,
        help="Project the cartesian velocity vector on the tangents of the internal"
             "coordinate. (in this case %s)" % name,
    )

def add_select_options(parser):
    parser.add_option(
        "-p", "--prefix", help="Format the indexes with the given prefix. The "
        "output will look like 'PREFIX.0000000 PREFIX.0000001 ...'"
    )
    parser.add_option(
        "--xyz", action='store_true', default=False,
        help="Append x, y and z extension to the prefixes. (Only applicable when "
        "-p is used. The output will look like 'PREFIX.0000000.x prefix.0000000.y "
        "PREFIX.0000000.z PREFIX.0000001.x ...'"
    )


