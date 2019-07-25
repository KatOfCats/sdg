Bootstrap: docker
From: gcc:6
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
cd /usr/local/
wget https://cmake.org/files/v3.4/cmake-3.4.1-Linux-x86_64.tar.gz
tar xvf cmake-3.4.1-Linux-x86_64.tar.gz
export PATH="`pwd`/cmake-3.4.1-Linux-x86_64/bin:$PATH"
rm -rf cmake-3.4.1-Linux-x86_64.tar.gz

# Install python3
apt-get install -y python3

#Install doxygen
apt-get install -y doxygen

#Install swig
apt-get install -y swig-4.0.0

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
