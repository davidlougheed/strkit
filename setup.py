#!/usr/bin/env python
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("./strkit/VERSION", "r") as vf:
    version = vf.read().strip()

setup(
    name="strkit",
    version=version,

    python_requires="~=3.9",
    install_requires=[
        "Flask>=2.2.5,<3.1",
        "orjson>=3.9.15,<4",
        "pysam>=0.19,<0.23",
        "numpy>=1.23.4,<1.27",
        "parasail>=1.2.4,<1.4",
        "scikit-learn>=1.2.1,<1.5",
        "scipy>=1.10,<1.14",
        "statsmodels>=0.14.0,<0.15",
        "strkit_rust_ext==0.12.1",
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

    packages=setuptools.find_namespace_packages(),
    include_package_data=True,

    entry_points={
        "console_scripts": ["strkit=strkit.entry:main"],
    },
)
