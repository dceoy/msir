#!/usr/bin/env python

from setuptools import setup, find_packages
from msir import __version__


setup(
    name='msir',
    version=__version__,
    description=(
        'Tandem repeat analyzer for microsatellite instability detection'
        ' by DNA-seq'
    ),
    packages=find_packages(),
    author='Daichi Narushima',
    author_email='dnarsil+github@gmail.com',
    url='https://github.com/dceoy/msir',
    include_package_data=True,
    install_requires=['biopython', 'docopt', 'pandas'],
    entry_points={'console_scripts': ['msir=msir.cli.main:main']},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Plugins',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    long_description="""\
msir
----

Tandem repeat analyzer for microsatellite instability detection by DNA-seq
"""
)
