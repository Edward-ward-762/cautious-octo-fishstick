#!/usr/bin/env nextflow

params.qc_input     = 'default'
params.mask_input   = 'default'
params.outdir       = './results'
params.publish_mode = 'copy'

process CAT_FASTA{
  tag "${meta.baseName}"
  label 'proces_low'

  publishDir "$params.outdir/$task.process", mode: "$params.publish_mode"

  input:
  tuple path(meta), path(ref_files)
  
  output:
  tuple val(meta), path("${meta.baseName}_cat.fa"), emit: cat_out

  script:
  """
  cat ${ref_files} > "${meta.baseName}_cat.fa"
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
    .unique { tuple -> [ tuple[0], tuple[1] ] }
    .groupTuple()

    CAT_FASTA(
      ch_mask.map{meta, ref -> [meta, ref] }
    )
    ch_cat_fasta = CAT_FASTA.out.cat_out

    ch_cat_fasta. map { fastq, cat_ref -> [ cat_ref, fastq ] }
      .view { cat_ref, fastq ->
        "Cat_Refs: $cat_ref, Fastq: $fastq"
      }

    ch_joined = ch_qc
      .join(ch_cat_fasta, by: [1])
      .meta {
        meta, fastq, cat ->
          if (cat) {
            [meta, fastq, cat]
          }
      }
      .view { meta, fastq, cat ->
        "Meta: $meta, Fastq: $fastq.baseName, Cat: $cat.baseName"
      }
}
