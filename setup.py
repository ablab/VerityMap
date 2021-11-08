import os
import sys
import subprocess


try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup
from distutils.command.build import build as DistutilsBuild
from distutils.spawn import find_executable

from veritymap.__version__ import __version__


# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)


description = \
    """
VerityMap aligns long PacBio Hi-Fi and ONT reads to
genome assemblies (including long repetitive regions)
"""


class MakeBuild(DistutilsBuild):
    def run(self):
        os.chdir(os.path.join(script_dir, "veritymap"))
        if not find_executable("make"):
            sys.exit("ERROR: 'make' command is unavailable")
        try:
            subprocess.check_call(["make"])
        except subprocess.CalledProcessError as e:
            sys.exit("Compilation error: ", e)
        os.chdir(script_dir)
        DistutilsBuild.run(self)


setup(
    name="VerityMap",
    version="2.0.0",
    description=description,
    url='https://github.com/ablab/VerityMap',
    author='Alla Mikheenko',
    author_email='al.miheenko@gmail.com',
    license='GNU General Public License v3.0',
    install_requires=[
    'plotly', 'python-slugify', 'biopython', 'numpy'],
    packages=['veritymap'],
    package_dir={'veritymap': 'veritymap'},
    package_data={'veritymap': ['build/bin/veritymap', 'config/*', '*', 'py_src/*', 'test_dataset/*',]},
    entry_points={
        'console_scripts': ['veritymap=veritymap.main:main']
    },
    cmdclass={'build': MakeBuild}
)
