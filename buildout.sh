#!/bin/bash
set -e

python=${1-"python"}

echo "### Blueprint RNAseq pipeline buildout  ###"
# set up binaries
echo "Download binaries..."
bintar=bp_pipeline_bin.tgz
wget -q http://genome.crg.es/~epalumbo/$bintar
tar xf $bintar
rm $bintar

# Install GEMtools
gemtools="GEMTools-static-i3-1.6.2"
if [ ! -d bin/gemtools-1.6.2-i3 ]; then
  echo "Install $gemtools"
  cd bin
  wget -q http://barnaserver.com/gemtools/releases/$gemtools.tar.gz
  tar xf $gemtools.tar.gz && rm $gemtools.tar.gz
  cd ..
fi

# Install Flux Capacitor
flux_ver="1.2.4"
flux="flux-capacitor-$flux_ver"
if [ ! -d bin/$flux ]; then
  echo "Install $flux"
  cd bin
  wget -q http://sammeth.net/artifactory/barna/barna/barna.capacitor/$flux_ver/$flux.tgz
  tar xf $flux.tgz && rm $flux.tgz
  cd ..
fi
exit 0
# set up virtual env
echo "Set up the Python virtualenv..."
virtualenv -q --no-site-packages -p $python env
echo "Activate the virtualenv..."
. env/bin/activate
echo "Install required Python packages..."
pip -q install -r pip-requirements.txt

# install RSeQC
echo "Download RSeQC 2.3.7..."
wget -q http://downloads.sourceforge.net/project/rseqc/RSeQC-2.3.7.tar.gz
tar xf RSeQC-2.3.7.tar.gz
rm RSeQC-2.3.7.tar.gz
cd RSeQC-2.3.7
echo "Install RSeQC 2.3.7..."
python setup.py install &> install.log
cd ..
cd env/lib/python2.7/site-packages/RSeQC-2.3.7-py2.7-linux-x86_64.egg
mv pysam pysam_old
ln -s ../pysam
rm -rf RSeQC-2.3.7
echo "Deactivate the virtualenv"
deactivate
echo "### done ###"
