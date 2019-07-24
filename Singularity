Bootstrap: docker
From: gcc:8
IncludeCmd: yes

%setup
# Copy any necessary files

%environment
 # Use bash as default shell
    SHELL=/bin/bash

    # Add paths
    PATH="/usr/local/bsg/bin:$PATH"

    # Export paths
    export PATH


%post
    # Make gpfs folder to hold mount
    mkdir /gpfs
    mkdir /shared
    mkdir /local
    mkdir /scratch
    # Create and move to build directory
    mkdir /root/build && cd /root/build

# run OS updates
apt-get update -y  && apt-get upgrade -y

# install cmake
apt-get install cmake

# Install python3
apt-get install python3

#Install doxygen
apt-get install doxygen

#Install swig
apt-get install swig-4.0.0

# Install sdg
\begin{python}
git clone https://github.com/bioinfologics/sdg
cd sdg
mkdir build
cd build
cmake -DBUILD_PYTHON_INTERFACE=ON ../
make
make install
\end{python}
