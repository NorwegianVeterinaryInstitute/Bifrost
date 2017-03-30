/*
* This code is part of the Bifrost package: https://github.com/karinlag/Bifrost
* This program should be run using nextflow. It makes heavy use
* of the program ariba for finding mlsts, antimicrobial resistance
* and virulence.
*/


Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set{ read_pairs1 }

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
    set pair_id1, file(reads) from read_pairs1
    file "mlst_db" from mlst_db

    output:
    file "${pair_id1}" into pair_id1

    """
    ariba run mlst_db/ref_db ${reads} ${pair_id1}

    """
}

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set{ read_pairs2 }


process run_ariba_amr_prep {
    publishDir params.out_dir + "/" + params.amr_results, mode: 'copy'

    output:
    file "out_amr_prepareref" into out_amr_prepareref

    """
    ariba getref ${params.amr_db} amr_db
    ariba prepareref -f amr_db.fa -m amr_db.tsv out_amr_prepareref
    """
}

process run_ariba_amr_pred {
    publishDir params.out_dir + "/" + params.amr_results, mode: 'copy'

    input:
    set pair_id2, file(reads) from read_pairs2
    file "out_amr_prepareref" from out_amr_prepareref

    output:
    file "${pair_id2}" into pair_id2

    """
    ariba run out_amr_prepareref ${reads} ${pair_id2}

    """
}