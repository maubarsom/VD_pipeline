#!/usr/bin/env nextflow

params.fastq_dir='/proj/virus/results/LCH/190402/preprocessing/'
params.ref_dir='/proj/virus/results/LCH/episome_seq/'

fq_files = Channel.fromFilePairs("${params.fastq_dir}/*/*.fq.gz",size: 3){ (it.name =~ /P[0-9]{3,5}_[0-9]{3,5}/)[0]}

ref_files = Channel.fromFilePairs("${params.ref_dir}/*.fasta", size: 1)

//fq_files.println()
//ref_files.println()

fq_files.into{bbmap_pe;
bbmap_unpe}

ref_files.into{bbmap_ref;
bbmap_unpe_ref}


process bbmap_paired{
  tag{"${sample_id} ${ref_id}"}
  publishDir "bbmap_results/${sample_id}", mode:'link'
   
  input:
  set sample_id,"reads_*.fq.gz" from bbmap_pe
  each ref from bbmap_ref

  output:
  set sample_id, ref_id, read_type, "${sample_id}_${ref_id}_${read_type}_basecov.txt" into python_data
  set sample_id, ref_id, read_type, "${sample_id}_${ref_id}_${read_type}.sam" into flagstat_data
  file "*.txt" into cov_files
  script:
  ref_id=ref[0]
  read_type="paired"
  """
  /proj/virus/tools/bbmap/38.44/bbmap.sh t="${task.cpus}" in1=reads_1.fq.gz in2=reads_2.fq.gz out="${sample_id}_${ref_id}_${read_type}.sam" covstats="${sample_id}_${ref_id}_${read_type}_covstat.txt" covhist="${sample_id}_${ref_id}_${read_type}_covhist.txt" basecov="${sample_id}_${ref_id}_${read_type}_basecov.txt" bincov="${sample_id}_${ref_id}_${read_type}_bincov.txt" ref="${ref[1][0]}"  
  """
}

process bbmap_unpaired{
  tag{"${sample_id} ${ref_id}"}
  publishDir "bbmap_results/${sample_id}", mode:'link'
   
  input:
  set sample_id,"reads_*.fq.gz" from bbmap_unpe
  each ref from bbmap_unpe_ref

  output:
  set sample_id, ref_id, read_type, "${sample_id}_${ref_id}_${read_type}_basecov.txt" into un_python_data
  set sample_id, ref_id, read_type, "${sample_id}_${ref_id}_${read_type}.sam" into un_flagstat_data
  file "*.txt" into un_cov_files
  script:
  ref_id=ref[0]
  read_type="unpaired"
  """
  /proj/virus/tools/bbmap/38.44/bbmap.sh t"=${task.cpus}" in=reads_3.fq.gz out="${sample_id}_${ref_id}_${read_type}.sam" covstats="${sample_id}_${ref_id}_${read_type}_covstat.txt" covhist="${sample_id}_${ref_id}_${read_type}_covhist.txt" basecov="${sample_id}_${ref_id}_${read_type}_basecov.txt" bincov="${sample_id}_${ref_id}_${read_type}_bincov.txt" ref="${ref[1][0]}"  
  """
}

mix_python_data = python_data.mix(un_python_data)

mix_flagstat_data = flagstat_data.mix(un_flagstat_data)

process flagstat{
  tag{"${sample_id} ${ref_id} ${read_type}"}
  publishDir "bbmap_results/${sample_id}", mode:'link'

  input:
  set sample_id, ref_id, read_type, bbmapOutput from mix_flagstat_data

  output:
  set sample_id, ref_id, read_type, "${sample_id}_${ref_id}_${read_type}_flagstat.txt" into flagstat_results

  script:
  """
   /proj/virus/tools/samtools/1.9/bin/samtools view -hbS "$bbmapOutput" | /proj/virus/tools/samtools/1.9/bin/samtools sort -O bam - | /proj/virus/tools/samtools/1.9/bin/samtools flagstat - > "${sample_id}_${ref_id}_${read_type}_flagstat.txt"
  """
}

process pythonPlot{
  tag{"${sample_id} ${ref_id} ${read_type}"}
  publishDir "bbmap_results/${sample_id}", mode:'link'

  input:
  set sample_id, ref_id, read_type, basecov from mix_python_data

  output:
  set sample_id, ref_id, read_type, "${sample_id}_${ref_id}_${read_type}_lineChartplot.png" into Plots

  script:
  """
  python3 /proj/virus/results/LCH/190405/gen_lineChart.py "${basecov}" "${sample_id}_${ref_id}_${read_type}_lineChartplot.png"
  """
}

process pythonGenTable{
tag{""}
publishDir "bbmap_results/", mode: 'link'

input:
file flagstat_files from flagstat_results.map{it[3]}.collect()

output:

file "LCH_episome_mapped_reads.tsv" into TableResult

script:
"""
python3 /proj/virus/results/LCH/190408/gen_table.py "./" "LCH_episome_mapped_reads.tsv"
"""

}





