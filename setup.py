from setuptools import setup, find_packages

PACKAGE = "inmembrane"
VERSION = "0.96.0-dev"  # keep in sync with inmembrane.__version__

DESCRIPTION = (
    "A modernized bioinformatics pipeline for proteome annotation and "
    "prioritization of bacterial surface / phage-receptor candidate proteins."
)
AUTHOR = "Andrew Perry, Bosco Ho, and contributors"
AUTHOR_EMAIL = "ajperry@pansapiens.com"
URL = "https://github.com/spyridon-ntokos/inmembrane"

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
    # The modern workflow relies on external CLI tools (SignalP 6.0, TMbed,
    # DeepLocPro, HMMER), which must be installed separately and configured
    # via inmembrane.config. We keep Python runtime dependencies minimal.
    install_requires=[
        # no mandatory non-stdlib dependencies for the core pipeline
        # (PyYAML, bs4, etc., can be installed separately if needed)
    ],
    license="BSD-3-Clause",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    zip_safe=False,
)
