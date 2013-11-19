#!/bin/bash
set -e

# set up binaries
bintar=bp_pipeline_bin.tgz
wget http://genome.crg.es/~epalumbo/$bintar
tar xf $bintar
rm $bintar

# set up virtual env
virtualenv --no-site-packages env
. env/bin/activate
pip install -r requirements.txt
deactivate

# install RSeQC
wget http://downloads.sourceforge.net/project/rseqc/RSeQC-2.3.7.tar.gz
tar xf RSeQC-2.3.7.tar.gz
rm RSeQC-2.3.7.tar.gz
cd RSeQC-2.3.7
python setup.py install
cd ..
