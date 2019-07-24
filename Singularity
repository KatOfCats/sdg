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
RUN apt-get install cmake

# Install python3
RUN apt-get install python3

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
