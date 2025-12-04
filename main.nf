#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def idr_output = params.idr_output ?: "idr_output"

process idr_call {
  tag "${pair}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${idr_output}", mode: 'copy'

  input:
    val pair
    path peaks1
    path peaks2

  script:
  """
  set -eux

  mkdir -p ${params.project_folder}/${idr_output}

  idr \\
    --samples ${peaks1} ${peaks2} \\
    --input-file-type narrowPeak \\
    --rank signal.value \\
    --output-file ${pair}_idr.txt \\
    --log-output-file ${pair}_idr.log \\
    --plot \\
    --idr-threshold ${params.idr_threshold}


  awk 'BEGIN{OFS="\\t"} !/^#/ {print \$1,\$2,\$3,"'${pair}'_IDR_peak_"NR,1000,".",\$7,\$8,\$9,\$10}' \\
      ${pair}_idr.txt > ${pair}_idr.narrowPeak
  """
}

workflow {

  def idr_csv = params.idr_pairs_csv ?: "idr_pairs.csv"

  rows = Channel
    .fromPath(idr_csv, checkIfExists: true)
    .splitCsv(header: true)

  rows = rows.filter { row ->
    ! file("${params.project_folder}/${idr_output}/${row.pair_name}_idr.narrowPeak").exists()
  }

  pairs_ch = rows.map { row ->
    tuple( row.pair_name, file(row.rep1_peaks), file(row.rep2_peaks) )
  }

  idr_call(pairs_ch)
}
