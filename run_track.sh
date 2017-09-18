#!/bin/bash -x

# Wrapper script for running nextflow pipelines more reproducibly

# Nextflow uses java, and we might run into memory issues. Thus
# we set the memory to something more managable

export NXF_OPTS='-Xms512M -Xmx2G'
script_directory=/home/karinlag/PycharmProjects/Bifrost

track_script=$1
template=$2
profile=$3
out_directory=$4

mkdir -p ${out_directory}/config_files
cp ${script_directory}/${track_script} ${out_directory}/config_files
cp ${template} ${out_directory}/config_files
cp ${script_directory}/conf/${profile}.config ${out_directory}/config_files

nextflow -c ${template} run -resume ${script_directory}/${track_script} -profile ${profile} --out_dir=${out_directory}
