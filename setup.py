
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

import subprocess, sys, os
from pathlib import Path

# https://stackoverflow.com/questions/20288711/post-install-script-with-python-setuptools
class DevelopCommand(develop):
    """Pre-installation for development mode."""
    def run(self):
        develop.run(self)

class InstallCommand(install):
    """Pre-installation for installation mode."""
    def run(self):
        install.run(self)

setup(name='PomatoData',
      version='0.1',
      description='PomatoData',
      author='Richard Weinhold',
      author_email='riw@wip.tu-berlin.de',
      url='https://github.com/richard-weinhold/pomato',
      packages=find_packages(),
      python_requires='>=3.6',
      include_package_data = True,
      install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'pathlib',
        'psutil',
        'geopandas',
        'geojson',
        'shapely',
        'requests',
        'scipy'
        ],
      cmdclass={
        'develop': DevelopCommand,
        'install': InstallCommand,}
     )

