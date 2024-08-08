FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y software-properties-common
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN add-apt-repository "deb http://security.ubuntu.com/ubuntu bionic-security main"
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    git \
    libeigen3-dev=3.3.4-4 \
    python3-dev=3.6.7-1~18.04 \
    python3-matplotlib=2.1.1-2ubuntu3 \
    gcc-7=7.5.0-3ubuntu1~18.04 \
    g++-7=7.5.0-3ubuntu1~18.04 \
    && rm -rf /var/lib/apt/lists/*

# Set up alternatives for gcc and g++
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-7 100 \
    && update-alternatives --set gcc /usr/bin/gcc-7 \
    && update-alternatives --set g++ /usr/bin/g++-7

ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++

# Verify compiler installation and paths
RUN gcc --version && g++ --version
RUN which gcc && which g++
RUN ls -l /usr/bin/gcc* /usr/bin/g++*

COPY postInstall /
RUN chmod +x /postInstall
RUN bash /postInstall

COPY cbf_tutorial /code/cbf_tutorial
COPY run /code/run_script
RUN chmod +x /code/run_script

RUN mkdir -p /results