#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, deicode development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools.command.egg_info import egg_info
from setuptools.command.develop import develop
from setuptools.command.install import install
import re
import ast
import os
from setuptools import find_packages, setup

# Dealing with Cython
USE_CYTHON = os.environ.get('USE_CYTHON', False)
ext = '.pyx' if USE_CYTHON else '.c'

# bootstrap numpy intall
# https://stackoverflow.com/questions/51546255/
# python-package-setup-setup-py-with-customisation
# -to-handle-wrapped-fortran


def custom_command():
    import sys
    if sys.platform in ['darwin', 'linux']:
        os.system('pip install numpy')


class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        custom_command()


class CustomDevelopCommand(develop):
    def run(self):
        develop.run(self)
        custom_command()


class CustomEggInfoCommand(egg_info):
    def run(self):
        egg_info.run(self)
        custom_command()


extensions = [
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Robust Aitchison compositional biplots from sparse count data')

with open('README.md') as f:
    long_description = f.read()

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('deicode/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

setup(name='deicode',
      version=version,
      license='BSD-3-Clause',
      description=description,
      long_description=long_description,
      author="deicode development team",
      author_email="cameronmartino@gmail.com",
      maintainer="deicode development team",
      maintainer_email="cameronmartino@gmail.com",
      packages=find_packages(),
      ext_modules=extensions,
      install_requires=[
          'numpy >= 1.12.1',
          'click',
          'pandas >= 0.10.0',
          'scipy >= 0.19.1',
          'nose >= 1.3.7',
          'scikit-bio > 0.5.3',
          'biom-format',
          'h5py', ],
      classifiers=classifiers,
      entry_points={
          'qiime2.plugins': ['q2-deicode=deicode.q2.plugin_setup:plugin'],
          'console_scripts':
              ['deicode=deicode.scripts._standalone_rpca:deicode']
      },
      # Inclusion of citations.bib in package_data based on how this is done in
      # q2-emperor's setup.py file
      package_data={'deicode': ['citations.bib']},
      cmdclass={'install': CustomInstallCommand,
                'develop': CustomDevelopCommand,
                'egg_info': CustomEggInfoCommand, },
      zip_safe=False)
