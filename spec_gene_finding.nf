/*
* This code is part of the Bifrost package: https://github.com/karinlag/Bifrost
* This program should be run using nextflow. It makes heavy use
* of the program ariba for finding mlsts, antimicrobial resistance
* and virulence.
*/


// First, define the input data that go into input channels
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into{ read_pairs_mlst; read_pairs_amr; read_pairs_vir }

// The following two processes are for MLST finding
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
    set pair_id, file(reads) from read_pairs_mlst
    file "mlst_db" from mlst_db

    output:
    file "${pair_id}" into pair_id_mlst

    """
    ariba run mlst_db/ref_db ${reads} ${pair_id}

    """
}

// These two processes are for AMR prediction
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
    set pair_id, file(reads) from read_pairs_amr
    file "out_amr_prepareref" from out_amr_prepareref

    output:
    file "${pair_id}" into pair_id_amr

    """
    ariba run out_amr_prepareref ${reads} ${pair_id}

    """
}


// These two processes are for virulence
process run_ariba_vir_prep {
    publishDir params.out_dir + "/" + params.vir_results, mode: 'copy'

    output:
    file "out_vir_prepareref" into out_vir_prepareref

    """
    ariba getref ${params.vir_db} vir_db
    ariba prepareref -f vir_db.fa -m vir_db.tsv out_vir_prepareref
    """
}

process run_ariba_vir_pred {
    publishDir params.out_dir + "/" + params.vir_results, mode: 'copy'

    input:
    set pair_id, file(reads) from read_pairs_vir
    file "out_vir_prepareref" from out_vir_prepareref

    output:
    file "${pair_id}" into pair_id_vir

    """
    ariba run out_vir_prepareref ${reads} ${pair_id}

    """
}