#!/usr/bin/env python3
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dryad",
    version="2.0.2",
    author="Kelsey Florek",
    author_email="kelsey.florek@slh.wisc.edu",
    description="A pipeline to construct reference free core-genome or SNP phylogenetic trees for examining prokaryote relatedness in outbreaks.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/k-florek/dryad",
    packages=setuptools.find_packages(),
    include_package_data=True,
    entry_points={
            "console_scripts":[
            'dryad = dryad_app.dryad:main',
            'dryad_report = dryad_app.dryad_report:main'
            ]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "pexpect>=4.8",
        ],
    python_requires='>=3.6',
)
