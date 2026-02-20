#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def idr_output        = params.idr_output        ?: "idr_output"
def pseudo_idr_output = params.pseudo_idr_output ?: "pseudo_idr_output"
def idr_rank          = params.idr_rank          ?: "signal.value"
def idr_input_type    = params.idr_input_type    ?: "narrowPeak"
def pseudo_peak_q     = params.qvalue            ?: 0.05

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

  idr \\
    --samples ${peaks1} ${peaks2} \\
    --input-file-type ${idr_input_type} \\
    --rank ${idr_rank} \\
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
  def macs_format = (params.seq == 'paired') ? 'BAMPE' : 'BAM'
  """
  set -eux

  samtools view -b -s 42.5 ${bam} -o ${rep_name}.pseudo1.bam
  samtools view -b -s 43.5 ${bam} -o ${rep_name}.pseudo2.bam

  samtools index ${rep_name}.pseudo1.bam
  samtools index ${rep_name}.pseudo2.bam

  # Call peaks on pseudo-replicates before IDR
  macs2 callpeak \\
    -t ${rep_name}.pseudo1.bam \\
    -n ${rep_name}_pseudo1 \\
    -f ${macs_format} \\
    -g ${params.macs3_genome} \\
    --outdir . \\
    --keep-dup all \\
    -q ${pseudo_peak_q}

  macs2 callpeak \\
    -t ${rep_name}.pseudo2.bam \\
    -n ${rep_name}_pseudo2 \\
    -f ${macs_format} \\
    -g ${params.macs3_genome} \\
    --outdir . \\
    --keep-dup all \\
    -q ${pseudo_peak_q}

  idr \\
    --samples ${rep_name}_pseudo1_peaks.narrowPeak ${rep_name}_pseudo2_peaks.narrowPeak \\
    --input-file-type ${idr_input_type} \\
    --rank ${idr_rank} \\
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
    path("${pair}_idr.sorted.chr.bed")
    path("${pair}_idr.sorted.narrowPeak")
    path("${pair}_idr.sorted.chr.narrowPeak")

  script:
  """
  set -euo pipefail

  cut -f1-3 ${idr_np} | sort -k1,1V -k2,2n > ${pair}_idr.sorted.bed

  awk 'BEGIN{OFS="\\t"} { if (\$1 !~ /^chr/) \$1 = "chr"\$1; print }' \
      ${pair}_idr.sorted.bed > ${pair}_idr.sorted.chr.bed

  sort -k1,1V -k2,2n ${idr_np} > ${pair}_idr.sorted.narrowPeak

  awk 'BEGIN{OFS="\\t"} { if (\$1 !~ /^chr/) \$1 = "chr"\$1; print }' \
      ${pair}_idr.sorted.narrowPeak > ${pair}_idr.sorted.chr.narrowPeak
  """
}


workflow {

  def idr_csv = params.idr_pairs_csv ?: "idr_pairs.csv"

  def rows = Channel
    .fromPath(idr_csv, checkIfExists: true)
    .splitCsv(header: true)

  def filtered_rows = rows.filter { row ->
    def bed_chr    = file("${params.project_folder}/${idr_output}/${row.pair_name}_idr.sorted.chr.bed")
    def narrow_chr = file("${params.project_folder}/${idr_output}/${row.pair_name}_idr.sorted.chr.narrowPeak")
    !(bed_chr.exists() && narrow_chr.exists())
  }

  def pairs_ch = filtered_rows.map { row ->
    assert row.pair_name && row.rep1_peaks && row.rep2_peaks : "idr_pairs.csv must contain: pair_name,rep1_peaks,rep2_peaks"
    def p1 = file(row.rep1_peaks.toString())
    def p2 = file(row.rep2_peaks.toString())
    assert p1.exists() : "rep1_peaks not found for ${row.pair_name}: ${p1}"
    assert p2.exists() : "rep2_peaks not found for ${row.pair_name}: ${p2}"
    tuple(row.pair_name.toString(), p1, p2)
  }

  def idr_results_ch = idr_call(pairs_ch)

  sort_idr_peaks(idr_results_ch)

  if( params.do_pseudo_idr ) {

    def pseudo_bam_csv = params.pseudo_idr_bam_csv ?: "pseudo_idr_bam.csv"

    def pseudo_rows = Channel
      .fromPath(pseudo_bam_csv, checkIfExists: true)
      .splitCsv(header: true)

    def pseudo_filtered_rows = pseudo_rows.filter { row ->
      ! file("${params.project_folder}/${pseudo_idr_output}/${row.rep_name}_pseudo_idr.narrowPeak").exists()
    }

    def pseudo_input_ch = pseudo_filtered_rows.map { row ->
      assert row.rep_name && row.bam : "pseudo_idr_bam.csv must contain: rep_name,bam"
      def b = file(row.bam.toString())
      assert b.exists() : "pseudo-IDR BAM not found for ${row.rep_name}: ${b}"
      tuple(row.rep_name.toString(), b)
    }

    pseudo_idr_from_bam(pseudo_input_ch)
  }
}
