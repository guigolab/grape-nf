#!/bin/bash
set -e

logfile="bp.buildout.log"
python=${1-"python"}

# set up binaries
echo "Downloading binaries..."
bintar=bp_pipeline_bin.tgz
wget -q http://genome.crg.es/~epalumbo/$bintar
tar xf $bintar
rm $bintar

# set up virtual env
echo "Setting up Python virtualenv..."
virtualenv -q --no-site-packages -p $python env
. env/bin/activate
pip -q install -r pip-requirements.txt
deactivate

# install RSeQC
echo "Installing RSeQC..."
wget i-q http://downloads.sourceforge.net/project/rseqc/RSeQC-2.3.7.tar.gz
tar xf RSeQC-2.3.7.tar.gz
rm RSeQC-2.3.7.tar.gz
cd RSeQC-2.3.7
python setup.py --quiet install
cd ..
echo "done"
