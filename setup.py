from setuptools import setup, find_packages

MAJOR = 0
MINOR = 1
PATCH = 0
VERSION = "{}.{}.{}".format(MAJOR, MINOR, PATCH)

#with open("pepr2ds/version.py", "w") as f:
#    f.write("__version__ = '{}'\n".format(VERSION))


setup(
    name='pepr2ds',
    version=VERSION,
    # url='http',
    license='GPL3',
    author='Thibault Tubiana and Phillippe Samer',
    author_email='T. Tubiana <tubiana.thibault@gmail.com>, P. Samer <samer@uib.no>',
    description='PePr2DS: Create and analysis dataset for peripheral membrane proteins',
    platforms=["Linux", "Solaris", "Mac OS-X", "darwin", "Unix", "win32"],
    # install_requires=['matplotlib',
    #                   'numpy',
    #                   'pandas',
    #                   'scipy >= 0.18',
    #                   'numba',
    #                   'mdtraj >= 1.9.5',
    #                   'PyQt5'],

    


    packages=find_packages(),
)