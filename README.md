# nf-idr

`nf-idr` is a Nextflow DSL2 module for reproducibility analysis of ChIP-seq peaks using IDR.

## What This Module Does

1. Runs regular IDR on replicate peak pairs from `idr_pairs.csv`.
2. Writes IDR outputs and sorted/chr-normalized peak files.
3. Optionally runs pseudo-IDR from BAM files when `do_pseudo_idr=true`.

## Input Files

### 1) Regular IDR pairs (`idr_pairs.csv`)

Required columns:

```text
pair_name,rep1_peaks,rep2_peaks
```

- `pair_name`: output prefix
- `rep1_peaks`, `rep2_peaks`: peak files (usually MACS narrowPeak)

### 2) Optional pseudo-IDR BAMs (`pseudo_idr_bam.csv`)

Used only when `--do_pseudo_idr true`.
Required columns:

```text
rep_name,bam
```

## Output

Under `${project_folder}/${idr_output}`:
- `${pair}_idr.txt`
- `${pair}_idr.log`
- `${pair}_idr.txt.png`
- `${pair}_idr.narrowPeak`
- `${pair}_idr.sorted.bed`
- `${pair}_idr.sorted.chr.bed`
- `${pair}_idr.sorted.narrowPeak`
- `${pair}_idr.sorted.chr.narrowPeak`

Pseudo-IDR outputs go to `${project_folder}/${pseudo_idr_output}`.

## Key Parameters

- `idr_pairs_csv`: path to regular IDR CSV
- `idr_threshold`: IDR cutoff (default: `0.05`)
- `idr_input_type`: IDR input type (default: `narrowPeak`)
- `idr_rank`: ranking column for IDR (default: `signal.value`)
- `do_pseudo_idr`: enable pseudo-IDR (default: `false`)
- `pseudo_idr_bam_csv`: path to pseudo-IDR BAM CSV
- `macs3_genome`: MACS genome size for pseudo peak calling (default: `mm`)
- `seq`: `paired` or `single` for pseudo peak calling

## Run

Regular IDR:

```bash
nextflow run main.nf -profile hpc --idr_pairs_csv idr_pairs.csv
```

Regular + pseudo-IDR:

```bash
nextflow run main.nf -profile hpc \
  --idr_pairs_csv idr_pairs.csv \
  --do_pseudo_idr true \
  --pseudo_idr_bam_csv pseudo_idr_bam.csv
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```
