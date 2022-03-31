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
    // Run Seurat pipeline


    emit:
    integration_ready = XXX.out
}

