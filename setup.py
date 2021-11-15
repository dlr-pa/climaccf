#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import find_packages, setup
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
  name='envlib',
  packages=find_packages(include=['envlib']),
  version='0.0.1',
  license='lgpl-3.0',
  description='Calculation of Algorithmic Climate Change Functions',
  #long_description=long_description,
  long_description_content_type='text/markdown',
  author='Abolfazl Simorgh',
  author_email='abolfazl.simorgh@uc3m.es',
  #url='',
  #download_url='',
  keywords=['Climate Impacts of Aviation', 'Algorithmic Climate Change Functions', 'Climate Hotspots', 'Non-CO2 ', 'Emissions'],
  install_requires=['xarray~=0.16.1', 'numpy==1.20.3', 'metpy~=1.1.0', 'h5netcdf~=0.8.1', 'scipy~=1.5.2', 'hyperopt~=0.2.5', 'pyparsing~=2.2.1', 'tomli~=1.0.2','setuptools~=49.6.0'],
  include_package_data=True,
  setup_requires=['pytest-runner==4.4'],
  tests_require=['pytest==4.4.1'],
  test_suite='tests',
  zip_safe=False,
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Education',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: Unix',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.8',
  ],
)
