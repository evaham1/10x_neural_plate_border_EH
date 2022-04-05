/*
 * Prepare RNA data for integration with ATAC
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.split                              = file("$baseDir/bin/seurat/split_seurat.R", checkIfExists: true)


params.split_options                                = [:]


// Include Seurat R processes
include {R as SPLIT} from "$baseDir/modules/local/r/main"                           addParams(  options:                        params.split_options,
                                                                                                script:                         analysis_scripts.split )


/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/

if(params.debug) {log.info Headers.build_debug_scripts_summary(analysis_scripts, params.monochrome_logs)}


/*------------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------------*/

workflow INTEGRATION_PREP {
    take:
    cell_cycle_data              //Channel: [[meta], [plot_dir, rds_dir]] - SEURAT_FILTERING.out.CELL_CYCLE.out
    transfer_labels              //Channel: [[meta], [plot_dir, rds_dir]]- TRANSFER_LABELS.out

    main:

    //
    // TRANSFER ORIGINAL CELL STATE LABELS:
    //

    // run a modified contamination_filt script to identify contaminating cell IDs and add them to metadata
    CONTAMINATION( cell_cycle_data )

    // use an edited transfer labels process to transfer labels of transfer_labels object into 'old' column of data
    ch_combined = CONTAMINATION.out
        .concat(transfer_labels)
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'tranfer_labels'], it] } // [[meta], [rds1, rds2, rds3, ...]]

    TRANSFER_OLD_LABELS( ch_combined )

    // Subset the input data to remove HH4
    SUBSET( TRANSFER_OLD_LABELS.out )
    CLUSTER_FULL( SUBSET.out )


    emit:
    integration_ready = XXX.out
}

