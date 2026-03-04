"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-chem-drugclip',  # Required
    version='4.1.4',  # Required
    description='Scipion plugin in order to use the DrugCLIP software',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-chem/scipion-chem-drugclip',  # Optional
    author='Blanca Pueche',  # Optional
    author_email='blanca.pueche@cnb.csic.es',  # Optional
    keywords='scipion drugclip scipion-3.0 cheminformatics',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={'pyworkflow.plugin': 'drugclip = drugclip'},
    package_data={  # Optional
       'drugclip': [ 'protocols.conf'],
    }
)
