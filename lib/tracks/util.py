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


__all__ = [
    "AtomFilter",
]


class AtomFilter(object):
    """A tool to test whether some atoms belong to a user defined set."""

    def __init__(self, filter_atoms=None):
        """Initialize the atom filter.

        The argument filter_atoms can be a list of atom indexes or a string with
        comma-separated atom indexes.
        """
        if isinstance(filter_atoms, str):
            self.filter_atoms = frozenset(int(word) for word in filter_atoms.split(","))
        elif filter_atoms is None:
            self.filter_atoms = None
        else:
            self.filter_atoms = frozenset(filter_atoms)

    def __call__(self, *test_indexes):
        """Test wither one of the indexes belongs to the predefined set."""
        if self.filter_atoms is None:
            return True
        return len(self.filter_atoms.intersection(test_indexes)) > 0

