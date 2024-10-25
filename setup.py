from setuptools import find_packages, setup
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='climaccf',
    packages=find_packages(include=['climaccf']),
    version='1.0',
    license='lgpl-3.0',
    description='Calculation of Algorithmic Climate Change Functions',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Abolfazl Simorgh',
    author_email='abolfazl.simorgh@uc3m.es',
    keywords=['Climate Impacts of Aviation', 'Algorithmic Climate Change Functions', 'Climate Hotspots', 'Non-CO2', 'Emissions'],
    install_requires=[
        'setuptools>=75.1.0', 
        'setuptools_scm~=8.1.0', 
        'pint~=0.24.3', 
        'xarray~=2024.10.0', 
        'numpy', 
        'metpy~=1.6.3', 
        'h5netcdf', 
        'scipy', 
        'hyperopt', 
        'pyparsing', 
        'pandas', 
        'casadi',  
        'tomli',
        'matplotlib', 
        'geojsoncontour', 
        'openpyxl', 
        'xlrd', 
        'PyYAML', 
        'pyproj',
        'pytest>=6.0' 
    ],
    include_package_data=True,
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
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)