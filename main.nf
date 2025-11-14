#!/usr/bin/env nextflow

params.qc_input   = 'default'
params.mask_input = 'default'

workflow {

  ch_qc = Channel.fromPath(params.qc_input)
    .splitCsv(header: true)
    .map { row ->
      [[id: row.sampled_id, genome_path: row.genome_path],row.fastq_path]
    }
    .view { meta, file ->
      "Sample: ${meta.id}, Genome: ${meta.genome_path}, Fastq: $file"
    }

  ch_mask = Channel.fromPath(params.mask_input)
    .splitCsv(header: true)
    .map { row ->
      [row.reads, row.roi_ref]
    }
    .view { reads, ref ->
      "Reads: $reads, Ref: $ref"
    }
}
