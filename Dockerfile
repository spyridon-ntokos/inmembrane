# Use an up-to-date Miniconda3 base (Debian 12/bookworm)
FROM continuumio/miniconda3:24.5.0-0

# Install system dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        tcsh \
        libxml2-dev \
        libxslt-dev \
        zlib1g-dev \
        python3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Configure conda channels
RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

# Install HMMER (the only truly required binary)
RUN conda install -y hmmer==3.1b2 && conda clean -a -y

# Set up working directory and environment
ENV APPS=/software/apps
WORKDIR ${APPS}
COPY . ${APPS}/inmembrane
WORKDIR ${APPS}/inmembrane

# Install Python dependencies and patch legacy deps
RUN pip install -U pip setuptools wheel && \
    # Fix old names in requirements.txt
    sed -i 's/BeautifulSoup/beautifulsoup4/' requirements.txt && \
    sed -i 's/^twill==0\.9\.1/twill==3.1.0/' requirements.txt && \
    pip install -r requirements.txt && \
    # Patch setup.py for modern libs
    sed -i "s/BeautifulSoup >= 3.2.1/beautifulsoup4>=4.11.1/" setup.py && \
    sed -i "s/BeautifulSoup/beautifulsoup4/g" setup.py && \
    sed -i "s/twill == 0.9.1/twill>=3.1.0/" setup.py && \
    pip install -e .

# Modernize the entrypoint script (print, exceptions, etc.)
RUN python -m lib2to3 -w /opt/conda/bin/inmembrane_scan

# Default entrypoint
# WORKDIR /data
ENTRYPOINT ["inmembrane_scan"]
