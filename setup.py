from setuptools import setup,find_packages
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name='femtoPy',
    version='0.1',
    description='Python code for data analysis and simulation of ultrafast phenomenon',
    url='http://github.com/dnp33/femtoPy',
    author='David N. Purschke',
    author_email='purschke@ualberta.ca',
    packages=['femtoPy','femtoPy/Diffusion','femtoPy/Bloch','femtoPy/THz','femtoPy/Optics','femtoPy/Relax'],
    keywords='femtosecond terahertz picosecond diffusion scattering-matrices',
)
