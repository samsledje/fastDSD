#!/usr/bin/env python

from setuptools import setup, find_packages
import fastdsd

setup(
    name="fastdsd",
    version=1,
    description="FastDSD",
    author="Enrico Maiorino",
    license="GPLv2",
    url="https://github.com/reemagit/DSD",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "fastdsd = fastdsd.__main__:main",
        ],
    },
    include_package_data=True,
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "networkx",
    ],
)
