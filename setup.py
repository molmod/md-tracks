#!/usr/bin/python
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

if __name__ == "__main__":
    from distutils.core import setup
    from glob import glob

    version = '0.001'

    setup(name='Tracks',
        version=version,
        description='Tracks provides tools for analyzing large trajectory files.',
        author='Toon Verstraelen',
        author_email='Toon.Verstraelen@UGent.be',
        url='https://molmod.ugent.be/zeobuilder',
        package_dir = {'': 'src'},
        py_modules = ['tracks'],
        scripts=glob("scripts/*"),
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python',
        ],
    )



