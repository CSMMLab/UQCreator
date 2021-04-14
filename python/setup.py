import os
import re
import sys
import platform
import subprocess
import multiprocessing

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

if sys.version_info < (3, 6):
    print("Python 3.6 or higher required, please upgrade.")
    sys.exit(1)

VERSION = "1.0"
REQUIREMENTS = ["numpy", "numba"]


setup(name='uqcreator',
      version=VERSION,
      author='Jonas Kusch, Jannick Wolters',
      description='UQCreator Python interface',
      long_description='',
      packages=["uqcreator",
                "uqcreator.closures",
                "uqcreator.problems"],
      install_requires=REQUIREMENTS,
      zip_safe=False)