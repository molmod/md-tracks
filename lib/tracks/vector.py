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


from tracks.core import load_track

import numpy


__all__ = [
    "TrackVector", "from_prefix",
    "dot", "cross", "triple",
    "bend", "dihed", "dist",
]


class TrackVector(object):
    def __init__(self, data, sub=slice(None)):
        self.data = data

    def __sub__(self, other):
        return TrackVector([c1 - c2 for c1, c2 in zip(self.data, other.data)])

    def __add__(self, other):
        return TrackVector([c1 - c2 for c1, c2 in zip(self.data, other.data)])

    def __mul__(self, other):
        if isinstance(other, TrackVector):
            return TrackVector([c1*c2 for c1, c2 in zip(self.data, other.data)])
        else:
            return TrackVector([c1*other for c1 in self.data])

    def __isub__(self, other):
        for c1, c2 in zip(self.data, other.data):
            c1 -= c2
        return self

    def __idiv__(self, other):
        if isinstance(other, TrackVector):
            for c1, c2 in zip(self.data, other.data):
                c1 /= c2
        else:
            for c1 in self.data:
                c1 /= other
        return self

    def norm(self):
        return numpy.sqrt(sum(c*c for c in self.data))


def from_prefix(prefix, sub=slice(None)):
    return TrackVector([load_track('%s.%s' % (prefix, c))[sub] for c in 'xyz'])


def dot(tv1, tv2):
    return sum(c1 * c2 for c1, c2 in zip(tv1.data, tv2.data))


def cross(tv1, tv2):
    return TrackVector([
        tv1.data[1] * tv2.data[2] - tv1.data[2] * tv2.data[1],
        tv1.data[2] * tv2.data[0] - tv1.data[0] * tv2.data[2],
        tv1.data[0] * tv2.data[1] - tv1.data[1] * tv2.data[0],
    ])


def triple(tv1, tv2, tv3):
    return (
        tv1.data[0] * (tv2.data[1] * tv3.data[2] - tv2.data[2] * tv3.data[1]) +
        tv1.data[1] * (tv2.data[2] * tv3.data[0] - tv2.data[0] * tv3.data[2]) +
        tv1.data[2] * (tv2.data[0] * tv3.data[1] - tv2.data[1] * tv3.data[0])
    )


def dist(p1, p2, v1=None, v2=None):
    """Compute the distance between two atoms at each time step."""
    pd = (p1 - p2).norm()
    if v1 is None:
        return pd
    else:
        vd = dot(v1-v2, p1-p2)/pd
        return pd, vd


def bend(v1, v2, v3, return_cos=False):
    """Compute the bending angle of three atoms at each time step."""
    delta_a = v1 - v2
    delta_b = v3 - v2
    # compute the norms
    norm_a = delta_a.norm()
    norm_b = delta_b.norm()
    # normalize the vectors
    delta_a /= norm_a
    delta_b /= norm_b
    # calculate the dot product
    cos = dot(delta_a, delta_b)
    cos = numpy.clip(cos, -1, 1)
    if return_cos: return cos
    angle = numpy.arccos(cos)
    return angle


def dihed(v1, v2, v3, v4, return_cos=False):
    """Compute the dihedral angle of four atoms at each time step."""
    delta_a = v1 - v2
    delta_b = v3 - v2
    delta_c = v4 - v3
    # compute the norm of b
    norm_b = delta_b.norm()
    # normalize the vector b
    delta_b /= norm_b
    # project a and c on the plane orthogonal to b
    dot_ab = dot(delta_a, delta_b)
    dot_cb = dot(delta_c, delta_b)
    tmp = delta_b*dot_ab
    delta_a -= delta_b*dot_ab
    delta_c -= delta_b*dot_cb
    # compute the norms of a' and c'
    norm_a = delta_a.norm()
    norm_c = delta_c.norm()
    # normalize the vectors a' and c'
    delta_a /= norm_a
    delta_c /= norm_c
    # calculate the dot product and the angle
    cos = dot(delta_a, delta_c)
    cos = numpy.clip(cos, -1, 1)
    if return_cos: return cos
    angle = numpy.arccos(cos)
    swap = (triple(delta_b, delta_a, delta_c) > 0)*2-1
    angle *= swap
    return angle


def oop(v1, v2, v3, v4):
    """Compute the distance from v4 to the plane defined by v1, v2 and v3 at each time step."""
    delta_a = v1 - v2
    delta_b = v3 - v2
    normal = cross(delta_a, delta_b)
    normal /= normal.norm()
    return dot(v4 - v2, normal)


def dtl(v1, v2, v3):
    """Compute the distance from v3 to the line defined by v1 and v2 at each time step."""
    delta_line = v1 - v2
    delta_line /= delta_line.norm()
    delta = v1 - v3
    return (delta - delta_line*dot(delta_line, delta)).norm()
