from setuptools import setup,find_packages

setup(
    name='femtoPy',
    version='0.1',
    description='Python code for data analysis and simulation of ultrafast phenomenon',
    url='http://github.com/dnp33/femtoPy',
    author='David N. Purschke',
    author_email='purschke@ualberta.ca',
    packages=['femtoPy','femtoPy/diffusion'],
    keywords='femtosecond terahertz thz picosecond diffusion',
    install_requires=['bandmat']
)
