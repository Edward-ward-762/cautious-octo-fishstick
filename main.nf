#!/usr/bin/env nextflow

params.qc_input   = 'default'
params.mask_input = 'default'

workflow {

  ch_qc = Channel.fromPath(params.qc_input)
    .splitCsv(header: true)
    .map { row ->
      [[id: row.sample_id, genome_path: row.genome_path],row.fastq_path]
    }
    .view { meta, file ->
      "Sample: ${meta.id}, Genome: ${meta.genome_path}, Fastq: $file"
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
}
