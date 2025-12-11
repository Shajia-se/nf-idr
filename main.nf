#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def idr_output        = params.idr_output        ?: "idr_output"
def pseudo_idr_output = params.pseudo_idr_output ?: "pseudo_idr_output"

process idr_call {
  tag "${pair}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${idr_output}", mode: 'copy'

  input:
    tuple val(pair), path(peaks1), path(peaks2)
  
  output:
    tuple val(pair), \
      path("${pair}_idr.narrowPeak"), \
      path("${pair}_idr.txt"), \
      path("${pair}_idr.log"), \
      path("${pair}_idr.txt.png")

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

process pseudo_idr_from_bam {
  tag "${rep_name}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${pseudo_idr_output}", mode: 'copy'

  input:
    tuple val(rep_name), path(bam)

  output:
    path("${rep_name}_pseudo_idr.narrowPeak")
    path("${rep_name}_pseudo_idr.txt")
    path("${rep_name}_pseudo_idr.log")
    path("${rep_name}_pseudo_idr.txt.png")

  script:
  """
  set -eux

  samtools view -b -s 42.5 ${bam} -o ${rep_name}.pseudo1.bam
  samtools view -b -s 43.5 ${bam} -o ${rep_name}.pseudo2.bam

  samtools index ${rep_name}.pseudo1.bam
  samtools index ${rep_name}.pseudo2.bam

  # GRCm39 effective genome size (deepTools)
  macs2 callpeak \\
    -t ${rep_name}.pseudo1.bam \\
    -n ${rep_name}_pseudo1 \\
    -f BAM \\
    -g 2654621783 \\
    --outdir . \\
    --keep-dup all \\
    -q ${params.qvalue}

  macs2 callpeak \\
    -t ${rep_name}.pseudo2.bam \\
    -n ${rep_name}_pseudo2 \\
    -f BAM \\
    -g 2654621783 \\
    --outdir . \\
    --keep-dup all \\
    -q ${params.qvalue}

  idr \\
    --samples ${rep_name}_pseudo1_peaks.narrowPeak ${rep_name}_pseudo2_peaks.narrowPeak \\
    --input-file-type narrowPeak \\
    --rank signal.value \\
    --output-file ${rep_name}_pseudo_idr.txt \\
    --log-output-file ${rep_name}_pseudo_idr.log \\
    --plot \\
    --idr-threshold ${params.idr_threshold}

  awk 'BEGIN{OFS="\\t"} !/^#/ {print \$1,\$2,\$3,"'${rep_name}'_PSEUDO_IDR_peak_"NR,1000,".",\$7,\$8,\$9,\$10}' \\
      ${rep_name}_pseudo_idr.txt > ${rep_name}_pseudo_idr.narrowPeak
  """
}

process sort_idr_peaks {

  tag "${pair}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${idr_output}", mode: 'copy'

  input:
    tuple val(pair), path(idr_np), path(idr_txt), path(idr_log), path(idr_png)

  output:
    path("${pair}_idr.sorted.bed")
    path("${pair}_idr.sorted.narrowPeak")

  script:
  """
  set -euo pipefail

  cut -f1-3 ${idr_np} | sort -k1,1V -k2,2n > ${pair}_idr.sorted.bed

  sort -k1,1V -k2,2n ${idr_np} > ${pair}_idr.sorted.narrowPeak
  """
}

workflow {

  def idr_csv = params.idr_pairs_csv ?: "idr_pairs.csv"

  rows = Channel
    .fromPath(idr_csv, checkIfExists: true)
    .splitCsv(header: true)

  rows = rows.filter { row ->
    def bed    = file("${params.project_folder}/${idr_output}/${row.pair_name}_idr.sorted.bed")
    def narrow = file("${params.project_folder}/${idr_output}/${row.pair_name}_idr.sorted.narrowPeak")
    !(bed.exists() && narrow.exists())
  }

  pairs_ch = rows.map { row ->
    tuple( row.pair_name, file(row.rep1_peaks), file(row.rep2_peaks) )
  }

  idr_results_ch = idr_call(pairs_ch)

  sort_idr_peaks(idr_results_ch)

  if( params.do_pseudo_idr ) {

    def pseudo_bam_csv = params.pseudo_idr_bam_csv ?: "pseudo_idr_bam.csv"

    pseudo_rows = Channel
      .fromPath(pseudo_bam_csv, checkIfExists: true)
      .splitCsv(header: true)

    pseudo_rows = pseudo_rows.filter { row ->
      ! file("${params.project_folder}/${pseudo_idr_output}/${row.rep_name}_pseudo_idr.narrowPeak").exists()
    }

    pseudo_input_ch = pseudo_rows.map { row ->
      tuple( row.rep_name, file(row.bam) )
    }

    pseudo_idr_from_bam(pseudo_input_ch)
  }
}
