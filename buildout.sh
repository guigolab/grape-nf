#!/bin/bash

log() {
    msg=$1
    date=`date`
    printf "$date - $msg\n"
}

set -e

python=${1-"python"}

echo "### Blueprint RNAseq pipeline buildout  ###"
# set up binaries
if [ ! -d bin ]; then
    log "Download and unpack binaries"
    bintar="bp_pipeline_bin_v1.22.tgz"
    wget -q http://genome.crg.es/~epalumbo/$bintar
    tar xf $bintar
    rm $bintar
fi

# Install GEMtools
gemtools="GEMTools-static-i3-1.6.2"
gemdir=`echo ${gemtools,,} | awk -F"-" '{print $1FS$4FS$3}'`
if [ ! -d bin/$gemdir ]; then
  log "Install $gemtools"
  cd bin
  wget -q http://barnaserver.com/gemtools/releases/$gemtools.tar.gz
  tar xf $gemtools.tar.gz && rm $gemtools.tar.gz
  cd ..
fi

# Install Flux Capacitor
flux_ver="1.2.4"
flux="flux-capacitor-$flux_ver"
if [ ! -d bin/$flux ]; then
  log "Install $flux"
  cd bin
  wget -q http://sammeth.net/artifactory/barna/barna/barna.capacitor/$flux_ver/$flux.tgz
  tar xf $flux.tgz && rm $flux.tgz
  cd ..
fi

# set up virtual env
env=${2-"."}
if [ ! -f $env/bin/python ]; then
    log "Set up the Python virtualenv"
    virtualenv -q --no-site-packages -p $python $env
    oldPath="PATH=\"\$VIRTUAL_ENV\/bin:\$PATH\""
    newPath="PATH=\"\$VIRTUAL_ENV\/bin:\$VIRTUAL_ENV\/bin\/$gemdir\/bin:\$VIRTUAL_ENV\/bin\/$flux\/bin:\$PATH\""
    sed -i 's/'$oldPath'/'$newPath'/g' $env/bin/activate
fi

# check if the virtual env is already active or activate it
if [[ ! $VIRTUAL_ENV ]];then
    log "Activate the virtualenv"
    . $env/bin/activate
fi

# install python packages
numpy=`python -c 'import numpy; print numpy.__version__' 2> /dev/null | cut -f1,2 -d.`
if [[ "`echo $numpy >= 1.7 | bc -l`" != 1 ]];then
    log "Install required Python packages"
    pip -q install -r pip-requirements.txt
fi

# download and unpack RSeQC
if [ ! -d RSeQC-2.3.7 ]; then
    log "Download and unpack RSeQC 2.3.7"
    wget -q http://downloads.sourceforge.net/project/rseqc/RSeQC-2.3.7.tar.gz
    tar xf RSeQC-2.3.7.tar.gz
    rm RSeQC-2.3.7.tar.gz
    cd RSeQC-2.3.7
fi

# instal RSeQC
if [ ! `python -c 'import bx; print bx.__version__' 2> /dev/null` ]; then
    log "Install RSeQC 2.3.7"
    python setup.py install &> install.log
    cd ..
    rm -rf RSeQC-2.3.7
    cd $env/lib/python2.7/site-packages/RSeQC-2.3.7-py2.7-linux-x86_64.egg
    mv pysam pysam_old
    ln -s ../pysam
fi

# deactivate the virtualenv if active
if [[ $VIRTUAL_ENV ]]; then
    log "Deactivate the virtualenv"
    deactivate
fi
log "### done ###"
