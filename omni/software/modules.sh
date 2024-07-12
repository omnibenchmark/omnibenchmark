#!/bin/bash
#
# Installs lua plus lmod (as module manager)
#
# Izaskun Mallona
# 27 May 2024


mkdir -p ~/soft/lua

cd $_

wget https://sourceforge.net/projects/lmod/files/lua-5.1.4.8.tar.gz/download

tar xzvf download

cd lua*8

./configure --with-static=yes --prefix=$HOME/soft/lua && make && make install

ldd ~/soft/lua/bin/lua

## only needed during lmod installation
export PATH=~/soft/lua/bin:$PATH

LMOD_VERS=8.7

mkdir -p ~/soft/lmod
cd $_

wget https://sourceforge.net/projects/lmod/files/Lmod-"$LMOD_VERS".tar.bz2/download

tar xfvj download && cd Lmod-"$LMOD_VERS"

./configure --prefix=$HOME/soft && make install

## also to .bashrc these three
export PATH=$HOME/soft/lmod/"$LMOD_VERS"/libexec:$PATH
source $HOME/soft/lmod/"$LMOD_VERS"/init/bash
export LMOD_CMD=$HOME/soft/lmod/"$LMOD_VERS"/libexec/lmod

lmod --version
