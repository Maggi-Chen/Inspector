FROM ubuntu:18.04


RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        bzip2 \
        curl \
        git \
        less \
        sudo \
        vim \
        wget \
        zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN sudo apt -y update
RUN sudo apt -y upgrade
RUN sudo apt -y install python2.7 python-pip
RUN pip install pysam
RUN python -m pip install -U matplotlib
RUN pip install statsmodels==0.10.1


RUN curl -L  https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 | tar -jxvf -
WORKDIR samtools-1.9
RUN ./configure --without-curses 
RUN make && make install
WORKDIR ..


RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.15/minimap2-2.15_x64-linux.tar.bz2 | tar -jxvf -
ENV PATH="minimap2-2.15_x64-linux/:${PATH}"

RUN pip install setuptools
RUN git clone https://github.com/fenderglass/Flye
WORKDIR Flye
RUN git checkout tags/2.8.3 -b inspector-flye
RUN python setup.py install
WORKDIR ..
RUN git clone https://github.com/Maggi-Chen/Inspector.git
ENV PATH="Inspector/:${PATH}"

