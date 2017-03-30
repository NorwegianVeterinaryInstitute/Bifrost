/*
* This code is part of the Bifrost package: https://github.com/karinlag/Bifrost
* This program should be run using nextflow. It makes heavy use
* of the program ariba for finding mlsts, antimicrobial resistance
* and virulence.
*/


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set{ read_pairs }

process run_ariba_mlst_prep {
    publishDir params.out_dir + "/" + params.mlst_results, mode: 'copy'

    output:
    file "mlst_db" into mlst_db

    """
    ariba pubmlstget "$params.mlst_scheme" mlst_db
    """
}

process run_ariba_mlst_pred {
    publishDir params.out_dir + "/" + params.mlst_results, mode: 'copy'

    input:
    set pair_id, file(reads) from read_pairs
    file "mlst_db" from mlst_db

    output:
    file "${pair_id}" into pair_id

    """
    ariba run mlst_db/ref_db ${reads} ${pair_id}

    """
}