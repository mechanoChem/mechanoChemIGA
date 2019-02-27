FROM ubuntu:16.04

#update ubuntu image with required packages
RUN apt-get update && apt-get -y install \
    cmake \
    g++ \
    gfortran \
    python \
    python-pip \
    openssh-client \
    openssh-server \
    liblapack-dev \
    libmpich-dev \
    libopenblas-dev \
    mpich \
    unzip \
    git \
    wget \
    emacs

#create directories
RUN  mkdir libraries && \
     chmod 777 \libraries

#run commands as non root
RUN adduser user
USER user

#install petsc
ENV PETSC_DIR=/libraries/petsc-3.8.3
ENV PETSC_ARCH=arch-linux2-c-debug
RUN cd /libraries && \
    wget -nc --quiet https://bitbucket.org/petsc/petsc/get/v3.8.3.tar.gz && \
    mkdir -p petsc-3.8.3 && tar -xf v3.8.3.tar.gz -C petsc-3.8.3 --strip-components 1 && \
    cd /libraries/petsc-3.8.3 && \
    ./configure --download-metis \
		--download-parmetis \
                --download-superlu_dist && \
    make all test && \
    rm /libraries/v3.8.3.tar.gz

#install PetIGA
RUN cd /libraries && \
    git clone https://bitbucket.org/dalcinl/PetIGA.git && \
    cd PetIGA && \
    git reset --hard 5ecf484 && \
    make all && \
    make test 
ENV PETIGA_DIR=/libraries/PetIGA

#install Sacado
RUN cd /libraries && \
    wget -nc --quiet http://trilinos.csbsju.edu/download/files/trilinos-11.8.1-Source.tar.gz && \
    tar -xzf trilinos-11.8.1-Source.tar.gz && \
    mkdir trilinos-11.8.1-Source/build && \
    mkdir trilinos && \
    cd /libraries/trilinos-11.8.1-Source/build && \
    cmake -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Teuchos=OFF -DCMAKE_INSTALL_PREFIX=/libraries/trilinos ../ && \
    make install && \
    cd /libraries && \
    rm -rf trilinos-11.8.1-Source && \
    rm /libraries/trilinos-11.8.1-Source.tar.gz
ENV TRILINOS_DIR=/libraries/trilinos

#download mechanoChem code
RUN cd /libraries && \
    git clone https://github.com/mechanoChem/mechanoChem.git

#install IGAkit (and dependancies) to convert .dat to .vtk
RUN pip install numpy && \
    pip install scipy && \
    pip install https://bitbucket.org/dalcinl/igakit/get/default.tar.gz