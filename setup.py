from setuptools import setup, find_packages

PACKAGE = "inmembrane"
DESCRIPTION = (
    "A modular bioinformatics pipeline for proteome annotation and "
    "prediction of bacterial surface-exposed proteins."
)
AUTHOR = "Andrew Perry, Bosco Ho, and contributors"
AUTHOR_EMAIL = "ajperry@pansapiens.com"
URL = "https://github.com/spyridon-ntokos/inmembrane-modern"
VERSION = "0.96.0"

setup(
    name=PACKAGE,
    version=VERSION,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    description=DESCRIPTION,
    packages=find_packages(
        include=[
            "inmembrane",
            "inmembrane.*",
        ]
    ),
    package_data={
        "inmembrane": [
            "protocols/*",
            "protocols/*/*",
        ]
    },
    scripts=["inmembrane_scan"],
    install_requires=[
        "beautifulsoup4>=4.11.1",
        "bs4",
        "cssselect",
        "lxml",
        "requests>=2.0.0",
        "semantic_version",
        "suds>=0.4",
        "twill>=3.1.0",
        "pyyaml",
    ],
    license="BSD-3-Clause",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    zip_safe=False,
)
