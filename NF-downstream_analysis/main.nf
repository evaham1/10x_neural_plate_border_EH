#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def skip_seurat_filtering = params.skip_seurat_filtering ? true : false
// def skip_scvelo = params.skip_scvelo ? true : false

/*-----------------------------------------------------------------------------------------------------------------------------
Log
-------------------------------------------------------------------------------------------------------------------------------*/
if(params.debug) {log.info Headers.build_debug_param_summary(params, params.monochrome_logs)}

/*------------------------------------------------------------------------------------*/
/* Workflow inclusions
--------------------------------------------------------------------------------------*/

include {METADATA} from "$baseDir/subworkflows/metadata/main"

include {SEURAT_FILTERING} from "$baseDir/subworkflows/seurat_filtering/main"                   addParams(  preprocessing_options:                  modules['preprocessing'],
                                                                                                            integration_options:                    modules['integration'],
                                                                                                            integration_qc_options:                 modules['integration_qc'],
                                                                                                            sex_filt_options:                       modules['sex_filt'],
                                                                                                            cell_cycle_options:                     modules['cell_cycle'],
                                                                                                            contamination_filt_options:             modules['contamination_filt'] )


// Modules and subworkflows for running scVelo/cellrank                                             
include {MERGE_LOOM} from "$baseDir/modules/local/merge_loom/main"                              addParams(  options:                                modules['merge_loom'] )

include {SEURAT_FILTERED_PROCESS} from "$baseDir/subworkflows/seurat_filtered_process/main"     addParams(  scatterplot3d_options:                  modules['scatterplot3d'],
                                                                                                            gene_module_options:                    modules['gene_modules'],
                                                                                                            state_classification_options:           modules['state_classification'],seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                            seurat_intersect_loom_options:          modules['seurat_intersect_loom'],
                                                                                                            scvelo_run_options:                     modules['scvelo_run'],
                                                                                                            gene_modules_latent_time_options:       modules['gene_modules_latent_time'])

// Subworkflows for subset stage, run and clusters from seurat object and run downstream analysis pipelines
include {SEURAT_STAGE_PROCESS} from "$baseDir/subworkflows/seurat_stage_process/main"           addParams(  split_options:                          modules['stage_split'],
                                                                                                            cluster_options:                        modules['stage_cluster'],
                                                                                                            gene_modules_options:                   modules['stage_gene_modules'],
                                                                                                            state_classification_options:           modules['stage_state_classification'],
                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                            seurat_intersect_loom_options:          modules['stage_seurat_intersect_loom'],
                                                                                                            scvelo_run_options:                     modules['stage_scvelo_run'],
                                                                                                            gene_modules_latent_time_options:       modules['stage_gene_modules_latent_time'])

include {SEURAT_RUN_PROCESS} from "$baseDir/subworkflows/seurat_run_process/main"               addParams(  split_options:                          modules['run_split'],
                                                                                                            cluster_options:                        modules['run_cluster'],
                                                                                                            gene_modules_options:                   modules['run_gene_modules'],
                                                                                                            state_classification_options:           modules['run_state_classification'],
                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                            seurat_intersect_loom_options:          modules['run_seurat_intersect_loom'],
                                                                                                            scvelo_run_options:                     modules['run_scvelo_run'],
                                                                                                            gene_modules_latent_time_options:       modules['run_gene_modules_latent_time'])

include {SEURAT_CLUSTERS_PROCESS} from "$baseDir/subworkflows/seurat_clusters_process/main"     addParams(  subset_options:                         modules['clusters_subset'],
                                                                                                            cluster_options:                        modules['clusters_cluster'],
                                                                                                            gene_modules_options:                   modules['clusters_gene_modules'],
                                                                                                            state_classification_options:           modules['clusters_state_classification'],
                                                                                                            seurat_h5ad_options:                    modules['seurat_h5ad'],
                                                                                                            seurat_intersect_loom_options:          modules['clusters_seurat_intersect_loom'],
                                                                                                            scvelo_run_options:                     modules['clusters_scvelo_run'],
                                                                                                            gene_modules_latent_time_options:       modules['clusters_gene_modules_latent_time'])


workflow {
    METADATA( params.input )

    /*------------------------------------------------------------------------------------*/
    /* Run inital seurat pipeline
    --------------------------------------------------------------------------------------*/
    // Set channel for cellranger counts
    METADATA.out
        .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
        .map {[it[0], it[1].collect{ file(it+"/cellranger/count/filtered_feature_bc_matrix", checkIfExists: true) }]}
        .set {ch_scRNAseq_counts}

    SEURAT_FILTERING( ch_scRNAseq_counts )
        
    /*------------------------------------------------------------------------------------*/
    /* Prepare inputs for scVelo
    --------------------------------------------------------------------------------------*/
   
    // Set channel for input looms
    METADATA.out
        .filter{ it[0].sample_id == 'NF-scRNAseq_alignment_out' }
        .map {[it[0], it[1].collect{ file(it+"/velocyto", checkIfExists: true) }]}
        .set {ch_loomInput}

    MERGE_LOOM( ch_loomInput )

    /*------------------------------------------------------------------------------------*/
    /* Run analysis on full filtered seurat object
    --------------------------------------------------------------------------------------*/
    SEURAT_FILTERED_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]})

    /*------------------------------------------------------------------------------------*/
    /* Split data and cluster batches and stages separately (+GMs)
    --------------------------------------------------------------------------------------*/ 
    SEURAT_STAGE_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]})  
    SEURAT_RUN_PROCESS( SEURAT_FILTERING.out.contamination_filt_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]})
    SEURAT_CLUSTERS_PROCESS( EXPLORATORY_ANALYSIS.out.state_classification_out, MERGE_LOOM.out.loom.map{it[1]}, SEURAT_FILTERING.out.annotations.map{it[1]})

}