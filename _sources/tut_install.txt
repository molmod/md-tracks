Installation
############


Disclaimer
==========

MD-Tracks is developed and tested in modern Linux environments. The
installation instructions below are given for a Linux system only. If you want
to use MD-Tracks on other operating systems such as Windows or OSX, you should
have a minimal computer geek status to get it working. We are always interested
in hearing from your installation adventures.


MolMod dependency
=================

`MolMod <http://molmod.github.com/molmod/>`_ is a Python library used by most
Python programs developed at the CMM. It must be installed before MD-Tracks can
be used or installed. Installation and download instructions can be found in the
`molmod documentation <http://molmod.github.com/molmod/tutorial/install.html>`_.
The instructions below only work if the MolMod package is installed.


External dependencies
=====================

Some software packages should be installed before MD-Tracks can be installed or
used. It is recommended to use the software package management of your Linux
distribution to install these dependencies.

The following software must be installed:

* Python 2.5, 2.6 or 2.7: http://www.python.org/
* Numpy >= 1.0: http://numpy.scipy.org/
* MatPlotLib >= 1.0: http://matplotlib.sourceforge.net/

Most Linux distributions can install this software with just a single terminal
command.

* Ubuntu 12.4::

    sudo apt-get install python python-numpy python-matplotlib

* Fedora 17::

    sudo yum install python numpy python-matplotlib


Installing the latest version of MD-Tracks
===========================================

The following series of commands will download the latest version of md-tracks,
and will then install it into your home directory. ::

    cd ~/build/
    git clone git://github.com/molmod/md-tracks.git
    (cd md-tracks; ./setup.py install --home=~)

You are now ready to start using MD-Tracks!


Upgrading to the latest version of MolMod and MD-Tracks
=======================================================

In case you want to upgrade MD-Tracks to the latests development version after
a previous install, then execute the following commands (in the same directory
that was originally used to install MD-Tracks)::

    cd ~/build/
    (cd molmod; git pull; rm -r ~/lib*/python/molmod*; ./setup.py install --home=~)
    (cd md-tracks; git pull; rm -r ~/lib*/python/tracks* ~/bin/tr-*; ./setup.py install --home=~)


Testing your installation
=========================

For the development and testing one needs to install additional packages:

* Nosetests >= 0.11: http://somethingaboutorange.com/mrl/projects/nose/0.11.2/
* Sphinx >= 1.0: http://sphinx.pocoo.org/

Most Linux distributions can install this software with just a single command:

* Ubuntu 12.4::

    sudo apt-get install python-nose python-sphinx

* Debian 5::

    su -
    apt-get install python-nose python-sphinx
    exit

* Fedora 17::

    sudo yum install python-nose sphinx

* Suse 11.2::

    sudo zypper install python-nose sphinx

Once these dependencies are installed, execute the following commands to run the
tests::

    cd ~/build/md-tracks
    nosetests -v

If some tests fail, post the output of the tests on the `MD-Tracks
mailing list <https://groups.google.com/forum/#!forum/MD-Tracks>`_.

