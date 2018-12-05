#!/bin/bash -x

# Wrapper script for running nextflow pipelines more reproducibly

# Nextflow uses java, and we might run into memory issues. Thus
# we set the memory to something more managable

export NXF_OPTS='-Xms512M -Xmx2G'
script_directory=$(dirname ${BASH_SOURCE[0]})

track_script=$1
template=$2
profile=$3
out_directory=$4

mkdir -p ${out_directory}/config_files
git --git-dir ${script_directory}/.git branch -v |grep "\*" | awk '{print $2, $3}' > ${out_directory}/config_files/pipeline_version.log
cp ${script_directory}/${track_script} ${out_directory}/config_files
cp ${template} ${out_directory}/config_files
cp ${script_directory}/conf/${profile}.config ${out_directory}/config_files

nextflow -c ${template} run -resume ${script_directory}/${track_script} -profile ${profile} --out_dir=${out_directory}
chmod -R 664 ${out_directory}
chmod -R a+X ${out_directory}
# spades sometimes gets weird permissions. Doing this to help removal
chmod -R 755 work 2> /dev/null
chmod -R -x+X work 2> /dev/null
