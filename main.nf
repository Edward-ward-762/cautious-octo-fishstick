#!/usr/bin/env nextflow

params.qc_input     = 'default'
params.mask_input   = 'default'
params.outdir       = './results'
params.publish_mode = 'copy'

process CAT_FASTA{
  tag "$meta"
  label 'proces_low'

  publishDir "$params.outdir/$task.process", mode: '$params.publish_mode'

  input:
  tuple val(meta), path(ref_files)
  
  output:
  tuple val(meta), path("${meta.baseName()}_cat.fa")

  script:
  """
  cat ${ref_files} > "${meta.baseName()}_cat.fa"
  """
}

workflow {

  ch_qc = Channel.fromPath(params.qc_input)
    .splitCsv(header: true)
    .map { row ->
      [[id: row.sample_id, genome_path: row.genome_path],row.fastq_path]
    }

  ch_mask = Channel.fromPath(params.mask_input)
    .splitCsv(header: true)
    .map { row ->
      [row.reads, row.roi_ref]
    }
    .groupTuple()

    ch_mask_2 = Channel.fromPath(params.mask_input)
    .splitCsv(header: true)
    .map { row ->
      [row.reads, row.roi_ref]
    }
    .unique { tuple -> [ tuple[0], tuple[1] ] }
    .groupTuple()
    .view { fastq, ref ->
      "Fastq: $fastq, Refs: $ref"
    }

    CAT_FASTA(
      ch_mask_2.map{meta, ref -> [meta, ref] }
    )
}
