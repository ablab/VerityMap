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

from tandemmapper2.__version__ import __version__


# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)


description = \
    """
TandemMapper2 aligns long PacBio Hi-Fi and ONT reads to
genome assemblies (including long repetitive regions)
"""


class MakeBuild(DistutilsBuild):
    def run(self):
        os.chdir(os.path.join(script_dir, "tandemmapper2"))
        if not find_executable("make"):
            sys.exit("ERROR: 'make' command is unavailable")
        try:
            subprocess.check_call(["make"])
        except subprocess.CalledProcessError as e:
            sys.exit("Compilation error: ", e)
        os.chdir(script_dir)
        DistutilsBuild.run(self)


setup(
    name="TandemMapper2",
    version="2.0.0",
    description=description,
    url='https://github.com/ablab/TandemMapper2',
    author='Alla Mikheenko',
    author_email='al.miheenko@gmail.com',
    license='GNU General Public License v3.0',
    install_requires=[
    'plotly', 'python-slugify', 'biopython', 'numpy'],
    packages=['tandemmapper2'],
    package_dir={'tandemmapper2': 'tandemmapper2'},
    package_data={'tandemmapper2': ['build/bin/tandem_mapper', 'config/*', '*', 'py_src/*', 'test_dataset/*',]},
    entry_points={
        'console_scripts': ['tandemmapper2=tandemmapper2.main:main']
    },
    cmdclass={'build': MakeBuild}
)
