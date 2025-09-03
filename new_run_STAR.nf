#!/usr/bin/env nextflow

/*
 * STAR paired-end alignment (single sample or many samples)
 * Usage examples:
 *   nextflow run new_.nf --fastq1 R1.fq.gz --fastq2 R2.fq.gz --sample_id SAMPLE --genomeDir /path/to/STAR/index
 *   nextflow run new_fixed_v3.nf --reads '/data/*_{R1,R2}.fastq.gz' --genomeDir /path/to/STAR/index
 */

nextflow.enable.dsl=2

// ----------------------------
// Parameters
// ----------------------------
params.fastq1    = null              // path to R1 (single-sample mode)
params.fastq2    = null              // path to R2 (single-sample mode)
params.sample_id = null              // sample name for single-sample mode
params.reads     = null              // glob for multi-sample mode: '/data/*_{R1,R2}.fastq.gz'
params.genomeDir = null              // REQUIRED: path to STAR genome index directory
params.outdir    = 'results2/star'    // where to publish outputs
params.threads   = 8                 // threads per alignment
params.rsemRef       = null     // REQUIRED to run RSEM 
params.forward_prob  = null     // optional: 0 (rev), 1 (fwd), 0.5 (unstranded). If null, RSEM estimates
params.rsem_extra    = null     // optional: any extra flags to pass to RSEM




// ----------------------------
// Sanity checks
// ----------------------------
if ( !params.genomeDir ) {
    log.error "Missing required parameter: --genomeDir </path/to/STAR/index>"
    System.exit(1)
}

if ( !params.reads && !(params.fastq1 && params.fastq2) ) {
    log.error "Provide either --reads '/glob/*_{R1,R2}.fastq.gz' for multiple samples OR --fastq1 and --fastq2 for a single sample."
    System.exit(1)
}

if ( params.reads && (params.fastq1 || params.fastq2) ) {
    log.warn "Both --reads and --fastq1/--fastq2 supplied. Proceeding with --reads and ignoring single-sample parameters."
}

// ----------------------------
// Build input channel
// ----------------------------
def read_pairs = ( params.reads )
    ? Channel.fromFilePairs( params.reads, flat: true, checkIfExists: true )
    : Channel.of( tuple( (params.sample_id ?: file(params.fastq1).simpleName) as String,
                         [ file(params.fastq1), file(params.fastq2) ] ) )

// ----------------------------
// Processes
// ----------------------------

process STAR_ALIGN {
  tag "${sample_id}"
  cpus params.threads
  publishDir "${params.outdir}/${sample_id}", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path("${sample_id}.Aligned.sortedByCoord.out.bam.bai"),                 emit: bai
    path("${sample_id}.Log.final.out"),                                     emit: log
    path("${sample_id}.ReadsPerGene.out.tab"),                              emit: counts
    tuple val(sample_id), path("${sample_id}.Aligned.toTranscriptome.out.bam"), emit: txbam  
    path("${sample_id}.*"),                                                 emit: all_files

  script:
    // define readCmd BEFORE the triple-quoted script
    def gz1 = reads[0].name.endsWith('.gz')
    def gz2 = reads[1].name.endsWith('.gz')
    def readCmd = (gz1 || gz2) ? '--readFilesCommand zcat' : ''

    """
    set -euo pipefail
    STAR \
      --genomeDir ${params.genomeDir} \
      --readFilesIn ${reads[0]} ${reads[1]} \
      ${readCmd} \
      --runThreadN ${task.cpus} \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --outFileNamePrefix ${sample_id}.

    samtools index -@ ${task.cpus} \
      -o ${sample_id}.Aligned.sortedByCoord.out.bam.bai \
         ${sample_id}.Aligned.sortedByCoord.out.bam
    """
}


// Sort the transcriptomic BAM
process SORT_TX_BAM_BY_NAME {
  tag "${sample_id}"
  cpus 2

  input:
    tuple val(sample_id), path(txbam)

  output:
    tuple val(sample_id), path("${sample_id}.tx.nameSorted.bam"), emit: txnamesort

  script:
  """
  set -euo pipefail
  samtools sort -@ ${task.cpus} -n -o ${sample_id}.tx.nameSorted.bam ${txbam}
  """
}



// Add after STAR_ALIGN process
process SAMTOOLS_INDEX {
  tag "${sample_id}"
  cpus 4
  input:
  tuple val(sample_id), path(bam)
  output:
  tuple val(sample_id), path("${bam}.bai")
  script:
  """
  set -euo pipefail
  samtools index -@ ${task.cpus} ${bam}
  """
}



//RSEM calculate
process RSEM_CALC {
  tag "${sample_id}"
  cpus params.threads
  publishDir "${params.outdir}/${sample_id}/rsem", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(txbam)

  output:
    path("${sample_id}.genes.results"),    emit: genes
    path("${sample_id}.isoforms.results"), emit: isoforms
    path("${sample_id}.*"),                emit: rsem_all

  script:
  """
  set -euo pipefail
  rsem-calculate-expression --bam\
    --no-bam-output \
    --paired-end \
    --num-threads ${task.cpus} \
    ${txbam} \
    ${params.rsemRef} \
    ${sample_id}
  """
}

// ----------------------------
// Workflow
// ----------------------------
workflow {
    main:
    STAR_ALIGN(read_pairs)
    RSEM_CALC( STAR_ALIGN.out.txbam )

    emit:
    bam    = STAR_ALIGN.out.bam
    logs   = STAR_ALIGN.out.log
    counts = STAR_ALIGN.out.counts
    files  = STAR_ALIGN.out.all_files
    bai    = STAR_ALIGN.out.bai
    txbam      = STAR_ALIGN.out.txbam
    rsem_genes = RSEM_CALC.out.genes
    rsem_iso   = RSEM_CALC.out.isoforms

}
