/*
 * Prepare RNA data for integration with ATAC
 */

/*------------------------------------------------------------------------------------*/
/* Set scripts used for downstream analysis
--------------------------------------------------------------------------------------*/

def analysis_scripts                                = [:]
analysis_scripts.contamination_ident = file("$baseDir/bin/seurat/6_contamination_filt.R", checkIfExists: true)
analysis_scripts.transfer_labels = file("$baseDir/bin/seurat/transfer_labels.R", checkIfExists: true)
analysis_scripts.subset                             = file("$baseDir/bin/seurat/subset_cells.R", checkIfExists: true)
analysis_scripts.cluster                            = file("$baseDir/bin/seurat/subset_cluster.R", checkIfExists: true)

params.contamination_ident_options   = [:]
params.transfer_labels_options   = [:]
params.subset_options = [:]
params.cluster_options                              = [:]

// Include Seurat R processes
include {R as CONTAMINATION_IDENT} from "$baseDir/modules/local/r/main"         addParams(      options: params.contamination_ident_options,
                                                                                                script: analysis_scripts.contamination_ident )
include {R as TRANSFER_LABELS_OLD} from "$baseDir/modules/local/r/main"         addParams(      options: params.transfer_labels_options,
                                                                                                script: analysis_scripts.transfer_labels )
include {R as SUBSET} from "$baseDir/modules/local/r/main"                      addParams(      options: params.subset_options,
                                                                                                script: analysis_scripts.subset )
include {R as CLUSTER} from "$baseDir/modules/local/r/main"                     addParams(      options: params.cluster_options,
                                                                                                script: analysis_scripts.cluster )


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

    // run contamination filt script with options to label rather than filter these cell states
    CONTAMINATION_IDENT( cell_cycle_data )

    // transfer labels process to transfer labels of transfer_labels object into 'old' column of data
    CONTAMINATION_IDENT.out
        .concat(transfer_labels)
        .map{it[1].findAll{it =~ /rds_files/}[0].listFiles()[0]}
        .collect()
        .map { [[sample_id:'tranfer_labels'], it] } // [[meta], [rds1, rds2, rds3, ...]]
        .set{ ch_combined }

    TRANSFER_LABELS_OLD( ch_combined )

    // Subset the input data to remove HH4
    SUBSET( TRANSFER_LABELS_OLD.out )
    CLUSTER_FULL( SUBSET.out )

    emit:
    integration_ready = CLUSTER_FULL.out
}

