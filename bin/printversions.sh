## This script attempts to print whatever versions we have of tools
## There are two scenarios, conda, and and not conda

profile=$1
output_file=$2


if [[ $profile == *"$conda"* ]];
    then
        conda list -n bifrost2022-fastqc >> ${output_file}
        conda list -n bifrost2022-multiqc >> ${output_file}
        conda list -n bifrost2022-bbtools >> ${output_file}
        conda list -n bifrost2022-trimmomatic >> ${output_file}
        conda list -n bifrost2022-spades >> ${output_file}
        conda list -n bifrost2022-bwa >> ${output_file}
        conda list -n bifrost2022-pilon >> ${output_file}
        conda list -n bifrost2022-prokka >> ${output_file}
        conda list -n bifrost2022-quast >> ${output_file}
        conda list -n bifrost2022-ariba >> ${output_file}

    else
        fastqc --version >> ${output_file}
        multiqc --version >> ${output_file}
        echo "bbduk.sh " `bbversion.sh` >> ${output_file}
        echo "trimmomatic " `trimmomatic -version` >> ${output_file}
        spades.py --version >> ${output_file}
        bwa &> tmpfile
        cat tmpfile |grep Version | awk '{print "bwa " $0}' >> ${output_file}
        rm tmpfile
        samtools --version |head  -1 >> ${output_file}
        pilon --version >> ${output_file} 
        prokka -v >> ${output_file}
        quast -v >> ${output_file}
        ariba version |head -1 >> ${output_file}


fi