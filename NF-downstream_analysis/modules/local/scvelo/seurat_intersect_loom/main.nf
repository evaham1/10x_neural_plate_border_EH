#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SEURAT_INTERSECT_LOOM {
    tag "$meta.sample_id"
    label 'process_medium'

    // maxForks 1

    
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[]) }

    container "alexthiery/10x-npb-scvelo:base-2.0.0"

    input:
        tuple val(meta), path(seurat), path(loom), path(annotations)

    output:
        tuple val(meta), path("*.loom"), emit: loom

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.prefix ? "${options.prefix}" : "${meta.sample_id}"
        
        """
        export HDF5_USE_FILE_LOCKING=FALSE
        $moduleDir/bin/seurat_intersect_loom.py --loomInput ${loom} --seuratInput ${seurat} --annotations ${annotations} --output ${prefix}_seurat_intersect.loom
        """
}