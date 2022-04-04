#!/usr/bin/env python
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="stronger",
    version="0.1.0-dev.10",

    python_requires="~=3.7",
    install_requires=[
        "pysam>=0.16.0.1,<0.20",
        "numpy>=1.20,<=1.22",
        "scikit-learn>=1.0,<1.1",
    ],

    description="A toolkit for analyzing variation in short(ish) tandem repeats.",
    long_description=long_description,
    long_description_content_type="text/markdown",

    url="https://github.com/davidlougheed/stronger",
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX",
    ],

    author="David Lougheed",
    author_email="david.lougheed@gmail.com",

    packages=setuptools.find_packages(),
    include_package_data=True,

    entry_points={
        "console_scripts": ["stronger=stronger.entry:main"],
    },
)
