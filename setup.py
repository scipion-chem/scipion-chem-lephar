"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import os, re

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

VERSION_FILE = os.path.join(here, "lephar", "version.py")

def get_version():
    with open(VERSION_FILE, "r") as f:
        content = f.read()
    match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', content, re.M)
    if not match:
        raise RuntimeError(f"Unable to find version string in {VERSION_FILE}")
    return match.group(1)
setup(
    name='scipion-chem-lephar',
    version=get_version(),
    description='Scipion-chem plugin for LePhar software use',
    long_description=long_description,
    url='https://github.com/scipion-chem/scipion-chem-lephar',
    author='Daniel Del Hoyo Gomez',
    author_email='ddelhoyo@cnb.csic.es',
    keywords='lephar ledock',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    package_data={
       'lephar': ['lephar_logo.jpg', 'protocols.conf'],
    },
    entry_points={
        'pyworkflow.plugin': 'lephar = lephar'
    }
)
