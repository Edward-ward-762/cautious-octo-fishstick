#!/usr/bin/env nextflow

params.qc_input     = 'default'
params.mask_input   = 'default'
params.outdir       = './results'
params.publish_mode = 'copy'

process CAT_FASTA{
  tag "${meta.id}"
  label 'proces_low'

  publishDir "$params.outdir/$task.process", mode: "$params.publish_mode"

  input:
  tuple val(meta), path(ref_files)
  
  output:
  tuple val(meta), path("${meta.id}_cat.fa"), emit: cat_out

  script:
  """
  cat ${ref_files} > "${meta.id}_cat.fa"
  """
}

workflow {

  //
  // CHANNEL: make channel from postplasmidsaurucQC samplesheet
  //

  ch_qc = Channel.fromPath(params.qc_input)
    .splitCsv(header: true)
    .map { row ->
      [[id: row.sample_id, genome_path: row.genome_path],row.fastq_path]
    }
    .view {meta, fastq ->
      "Meta: $meta, Fastq: $fastq"
    }

  //
  // CHANNEL: make channel from cas9point4 master samplesheet
  //

  ch_mask = Channel.fromPath(params.mask_input)
    .splitCsv(header: true)
    .map { row ->
      [row.reads, row.roi_ref]
    }
    .unique { tuple -> [ tuple[0], tuple[1] ] }
    .groupTuple()
  
  ch_mask = ch_mask
    .map{ fastq, ref ->
      [ref, fastq]
    }

  ch_join = ch_qc
    .join(ch_mask, by: [1])
    .map{ fastq, meta, ref -> 
      meta.roi_ref = ref
      [meta, fastq]
    }
    .view()

  //
  // MODULE: cat references fasta files together
  //

  CAT_FASTA(
    ch_join.map{meta, ref -> [meta, meta.roi_ref] }
  )
  ch_cat_fasta = CAT_FASTA.out.cat_out

}
