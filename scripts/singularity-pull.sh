#!/bin/bash
set -e
set -u

imgUrls=(
    grapenf/ihec:latest
)

for image in ${imgUrls}; do
    name=${image//[:\/]/-}.img
    singularity pull --name ${name} docker://$img
done