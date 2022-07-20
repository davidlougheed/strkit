#!/usr/bin/env python
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="strkit",
    version="0.3.0.rc1",

    python_requires="~=3.7",
    install_requires=[
        "Flask>=2.1.3,<2.2",
        "pysam>=0.16.0.1,<0.20",
        "numpy>=1.20,<=1.22",
        "parasail>=1.2.4,<1.3",
        "scikit-learn>=1.0,<1.2",
    ],

    description="A toolkit for analyzing variation in short(ish) tandem repeats.",
    long_description=long_description,
    long_description_content_type="text/markdown",

    url="https://github.com/davidlougheed/strkit",
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
        "console_scripts": ["strkit=strkit.entry:main"],
    },
)
