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

    def __div__(self, other):
        if isinstance(other, TrackVector):
            return TrackVector([c1/c2 for c1, c2 in zip(self.data, other.data)])
        else:
            return TrackVector([c1/other for c1 in self.data])

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


def bend(p1, p2, p3, v1=None, v2=None, v3=None, return_cos=False):
    """Compute the bending angle of three atoms at each time step."""
    p_delta_a = p1 - p2
    p_delta_b = p3 - p2
    # compute the norms
    p_norm_a = p_delta_a.norm()
    p_norm_b = p_delta_b.norm()
    # calculate the cosine
    p_cos = dot(p_delta_a, p_delta_b)/p_norm_a/p_norm_b
    p_cos = numpy.clip(p_cos, -1, 1)
    if v1 is None:
        if return_cos: return p_cos
        p_angle = numpy.arccos(p_cos)
        return p_angle
    else:
        v_norm_a = dot(v1-v2, p_delta_a)/p_norm_a
        v_norm_b = dot(v3-v2, p_delta_b)/p_norm_b

        t = dot(p_delta_a, p_delta_b)
        n = p_norm_a*p_norm_b
        dt = dot(v1, p_delta_b) + dot(v3, p_delta_a) - dot(v2, p_delta_a) - dot(v2, p_delta_b)
        dn = v_norm_a*p_norm_b + v_norm_b*p_norm_a

        v_cos = (n*dt-t*dn)/n**2
        if return_cos: return p_cos, v_cos
        p_angle = numpy.arccos(p_cos)
        v_angle = -1/numpy.sqrt(1 - p_cos**2)*v_cos
        return p_angle, v_angle



def dihed(p1, p2, p3, p4, v1=None, v2=None, v3=None, v4=None, return_cos=False):
    """Compute the dihedral angle of four atoms at each time step."""
    p_delta_a = p1 - p2
    p_delta_b = p3 - p2
    p_delta_c = p4 - p3
    # compute the norm of b
    p_norm_b = p_delta_b.norm()
    # normalize the vector b
    p_normed_b = p_delta_b/p_norm_b
    # project a and c on the plane orthogonal to b
    p_dot_ab = dot(p_delta_a, p_normed_b)
    p_dot_cb = dot(p_delta_c, p_normed_b)
    p_proj_a = p_delta_a - p_normed_b*p_dot_ab
    p_proj_c = p_delta_c - p_normed_b*p_dot_cb
    # compute the norms of a' and c'
    p_proj_norm_a = p_proj_a.norm()
    p_proj_norm_c = p_proj_c.norm()
    # calculate the cosine and/or the angle
    p_t = dot(p_proj_a, p_proj_c)
    p_n = p_proj_norm_a*p_proj_norm_c
    p_cos = p_t/p_n
    p_cos = numpy.clip(p_cos, -1, 1)
    if not return_cos:
        p_angle = numpy.arccos(p_cos)
        swap = (triple(p_delta_b, p_delta_a, p_delta_c) > 0)*2-1
        p_angle *= swap
    if v1 is None:
        if return_cos:
            return p_cos
        else:
            return p_angle
    else:
        v_delta_a = v1 - v2
        v_delta_b = v3 - v2
        v_delta_c = v4 - v3
        v_norm_b = dot(p3-p2,v3-v2)/p_norm_b
        v_normed_b = (v_delta_b - p_normed_b*v_norm_b)/p_norm_b
        v_proj_a = v_delta_a - v_normed_b*dot(p_normed_b,p_delta_a) \
                             - p_normed_b*dot(v_normed_b,p_delta_a) \
                             - p_normed_b*dot(p_normed_b,v_delta_a)
        v_proj_c = v_delta_c - v_normed_b*dot(p_normed_b,p_delta_c) \
                             - p_normed_b*dot(v_normed_b,p_delta_c) \
                             - p_normed_b*dot(p_normed_b,v_delta_c)
        v_proj_norm_a = dot(p_proj_a,v_proj_a)/p_proj_norm_a
        v_proj_norm_c = dot(p_proj_c,v_proj_c)/p_proj_norm_c
        v_t = dot(v_proj_a, p_proj_c) + dot(p_proj_a, v_proj_c)
        v_n = p_proj_norm_a*v_proj_norm_c + v_proj_norm_a*p_proj_norm_c
        v_cos = (p_n*v_t - v_n*p_t)/p_n**2
        if return_cos:
            return p_cos, v_cos
        else:
            v_angle = -1/numpy.sqrt(1 - p_cos**2)*v_cos*swap
            return p_angle, v_angle


def oop(p1, p2, p3, p4, v1=None, v2=None, v3=None, v4=None):
    """Compute the distance from p4 to the plane defined by p1, p2 and p3 at each time step."""
    p_delta_a = p1 - p2
    p_delta_b = p3 - p2
    p_ortho_1 = p4 - p1

    p_norm_a = p_delta_a.norm()
    p_normed_a = p_delta_a/p_norm_a
    p_ortho_2 = p_ortho_1 - p_normed_a*dot(p_normed_a,p_ortho_1)

    p_proj_b = p_delta_b - p_normed_a*dot(p_normed_a,p_delta_b)
    p_norm_proj_b = p_proj_b.norm()
    p_normed_proj_b = p_proj_b/p_norm_proj_b
    p_ortho_3 = p_ortho_2 - p_normed_proj_b*dot(p_normed_proj_b,p_ortho_2)
    p_oop = p_ortho_3.norm()

    swap = (triple(p_delta_a, p_delta_b, p_ortho_1) > 0)*2-1
    p_oop *= swap
    if v1 is None:
        return p_oop
    else:
        v_delta_a = v1 - v2
        v_delta_b = v3 - v2
        v_ortho_1 = v4 - v1

        v_norm_a = dot(p_delta_a, v_delta_a)/p_norm_a
        v_normed_a = (v_delta_a - p_normed_a*v_norm_a)/p_norm_a
        v_ortho_2 = v_ortho_1 - v_normed_a*dot(p_normed_a,p_ortho_1) \
                              - p_normed_a*dot(v_normed_a,p_ortho_1) \
                              - p_normed_a*dot(p_normed_a,v_ortho_1)

        v_proj_b = v_delta_b - v_normed_a*dot(p_normed_a,p_delta_b) \
                             - p_normed_a*dot(v_normed_a,p_delta_b) \
                             - p_normed_a*dot(p_normed_a,v_delta_b)
        v_norm_proj_b = dot(p_proj_b, v_proj_b)/p_norm_proj_b
        v_normed_proj_b = (v_proj_b - p_normed_proj_b*v_norm_proj_b)/p_norm_proj_b
        v_ortho_3 = v_ortho_2 - v_normed_proj_b*dot(p_normed_proj_b,p_ortho_2) \
                              - p_normed_proj_b*dot(v_normed_proj_b,p_ortho_2) \
                              - p_normed_proj_b*dot(p_normed_proj_b,v_ortho_2)
        v_oop = dot(p_ortho_3, v_ortho_3)/p_oop

        return p_oop, v_oop


def dtl(v1, v2, v3):
    """Compute the distance from v3 to the line defined by v1 and v2 at each time step."""
    delta_line = v1 - v2
    delta_line /= delta_line.norm()
    delta = v1 - v3
    return (delta - delta_line*dot(delta_line, delta)).norm()
