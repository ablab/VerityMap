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

from tandemmapper.__version__ import __version__


# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)


requirements_fn = os.path.join(script_dir, 'requirements.txt')
requirements = []
with open(requirements_fn) as f:
    for line in f:
        line = line.strip()
        requirements.append(line)


description = \
    """
TandemMapper2 aligns long PacBio Hi-Fi and ONT reads to
genome assemblies (including long repetitive regions)
"""


class MakeBuild(DistutilsBuild):
    def run(self):
        if not find_executable("make"):
            sys.exit("ERROR: 'make' command is unavailable")
        try:
            subprocess.check_call(["make"])
        except subprocess.CalledProcessError as e:
            sys.exit("Compilation error: ", e)
        DistutilsBuild.run(self)


setup(
    name="TandemMapper2",
    version=__version__,
    description=description,
    url='https://github.com/ablab/TandemMapper2',
    author='Alla Mikheenko',
    author_email='al.miheenko@gmail.com',
    license='GNU General Public License v3.0',
    install_requires=requirements,
    packages=['tandemmapper2'],
    package_dir={'tandemmapper2': 'tandemmapper2'},
    package_data={'tandemmapper2': ['build/bin/tandem_mapper', 'config/*', '*', 'py_src/*', 'test_dataset/*',]},
    entry_points={
        'console_scripts': ['tandemmapper2=tandemmapper2.tandemmapper2:main']
    },
    cmdclass={'build': MakeBuild}
)
