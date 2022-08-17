#!/bin/bash

# Wrapper script for running nextflow pipelines more reproducibly

# Nextflow uses java, and we might run into memory issues. Thus
# we set the memory to something more managable

export NXF_OPTS='-Xms512M -Xmx2G'
script_directory=$(dirname ${BASH_SOURCE[0]})

track_script=$1
template=$2
profile=$3
out_directory=$4
workdir=${5:-$USERWORK/bifrost_work}

#to add to report name
now=$(date +"%Y%m%d_%H%M")

mkdir -p ${out_directory}/config_files
git --git-dir ${script_directory}/.git branch -v |grep "\*" | awk '{print $2, $3}' > ${out_directory}/config_files/pipeline_version.log
bash ${script_directory}/bin/printversions.sh ${profile} ${out_directory}/config_files/software_versions.txt
cp ${script_directory}/${track_script} ${out_directory}/config_files
cp ${template} ${out_directory}/config_files


echo "TEMPORARY WORKING DIRECTORY IS ${workdir}"
nextflow -c ${template} run -resume ${script_directory}/${track_script} \
-profile ${profile} --out_dir=${out_directory} -work-dir ${workdir} -with-report ${now}_run_report.html
