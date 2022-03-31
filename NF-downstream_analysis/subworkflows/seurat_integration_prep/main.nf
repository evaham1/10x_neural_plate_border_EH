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
    new_labeled_data             //Channel: [[meta], [plot_dir, rds_dir]]- SEURAT_STAGE_PROCESS_CONTAM.out

    main:

    //
    // TRANSFER OVER OLD CELL STATE LABELS:
    //

    // run a modified contamination_filt script to identify contaminating cell IDs and add them to metadata
    CONTAMINATION( cell_cycle_data )

    // use an edited transfer labels process to transfer labels of transfer_labels object into 'old' column of data
    //channel operation similar to the one below to read in both CONTAMINATION.out and TRANSFER LABELS object
    //TRANSFER_OLD_LABELS(transfer_labels)

    // Subset the input data to remove HH4
    SUBSET( TRANSFER_OLD_LABELS.out )
    CLUSTER_FULL( SUBSET.out )

    //
    // COMBINE OLD AND NEW LABELS:
    //

    // Transfer labels from individual stages to merged data
    ch_labels = STATE_CLASSIFICATION.out
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'all_stages'], it] } // [[meta], [rds1, rds2, rds3, ...]]
        .combine( CLUSTER_FULL.out ) //[[sample_id:all_stages], [HH7, ss8, HH6, ss4, HH4, HH5], [sample_id:NF-scRNA-input], [rds_files, plots]]
        .map{[it[0], it[1] + it[3]]} //[[sample_id:all_stages], [HH7, ss8, HH6, ss4, HH4, HH5, rds_files, plots]
        //.view() //[[sample_id:all_stages], [HH6, HH4, ss8, ss4, HH7, HH5, cell_cycle_data.RDS]]
    TRANSFER_LABELS( ch_labels )


    emit:
    integration_ready = XXX.out
}

