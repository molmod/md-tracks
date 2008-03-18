# Tracks provides tools for analyzing large trajectory files.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
