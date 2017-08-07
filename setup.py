"""
materials-toolbox: making HPC a little easier
"""

from os.path import abspath, dirname
from setuptools import setup, find_packages

project_dir = abspath(dirname(__file__))

setup(
    name='materials-toolbox',
    version='1.0.0',
    description='Convenience tools for working with comp chem HPC',
    long_description="""
Provides some essential tools for working with computational chemistry HPC.
""",
    url="https://github.com/utf/materials-toolbox",
    author="Alex Ganose",
    author_email="alexganose@googlemail.com",
    license='MIT',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Awk',
        'Programming Language :: Unix Shell',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords='chemistry ase dft vasp',
    packages=find_packages(),
    install_requires=['ase', 'spglib', 'numpy', 'pymatgen', 'seekpath'],
    entry_points={'console_scripts': [
                      'sym = matools.sym:main',
                      'conv = matools.conv:main',
                      'prim = matools.prim:main',
                      'super = matools.super:main',
                      'mp-get = matools.mp:main',
                      'reorient-poscar = matools.reorient:main']},
    scripts=['bin/cpos',
             'bin/bandgap']
    )
