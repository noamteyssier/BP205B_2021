#!/usr/bin/env bash

# initialize DeMultiSeq submodule
git submodule init
git submodule update

cd src/DeMultiSeq
./setup.sh
cd -
