/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMicroc.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Modify fasta channel to include meta data
if (params.fasta) { ch_fasta =  Channel.fromPath(params.fasta) } else { exit 1, 'Fasta reference genome not specified!' }
ch_fasta_meta = ch_fasta.map{ it -> [[id:it[0].baseName], it] }.collect()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
//include { FASTQ_ALIGN_BWA } from '../subworkflows/nf-core/fastq_align_bwa/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'
//include { PAIRTOOLS_PARSE             } from '../modules/nf-core/pairtools/parse/main'
//include { PAIRTOOLS_SORT              } from '../modules/nf-core/pairtools/sort/main'
//include { PAIRTOOLS_DEDUP             } from '../modules/nf-core/pairtools/dedup/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { BWA_MEM                     } from '../modules/nf-core/bwa/mem/main'
include { BWA_ALN                     } from '../modules/nf-core/bwa/aln/main'
include { BWA_SAMSE                   } from '../modules/nf-core/bwa/samse/main'

include { PAIRIX                      } from '../modules/nf-core/pairix/main'
include { COOLER_CLOAD                } from '../modules/nf-core/cooler/cload/main'
include { HICEXPLORER_HICPCA          } from '../modules/nf-core/hicexplorer/hicpca/main'


include { PAIRTOOLS_PARSE             } from '../modules/local/pairtools/parse/main'
include { PAIRTOOLS_SORT              } from '../modules/local/pairtools/sort/main'
include { PAIRTOOLS_DEDUP             } from '../modules/local/pairtools/dedup/main'
include { PAIRTOOLS_SPLIT             } from '../modules/local/pairtools/split/main'
include { PAIRTOOLS_QC                } from '../modules/local/pairtools/pairtools_qc/main'

include { CHROMSIZES                  } from '../modules/local/chromsizes'
include { JUICER_PRE                  } from '../modules/local/juicer/pre.nf'
include { JUICER_PRE     as   JUICER_PRE_250k  } from '../modules/local/juicer/pre.nf'
include { JUICER_PRE     as   JUICER_PRE_50k  } from '../modules/local/juicer/pre.nf'
//SUB-WORKFLOWS
//include { CAT_SAMPLES                 } from '../subworkflows/local/cat_samples.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MICROC {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    SAMTOOLS_FAIDX (
            ch_fasta_meta,
            [[], []]
    )

    //INPUT_CHECK.out.reads.view()

    ch_nocat = INPUT_CHECK.out.reads
                     .map{ meta, path -> 
                        id=meta.subMap('id')
                        meta=meta
                        path=path
                        [id.id, meta, path]
                      }
                     .groupTuple()
                     .filter{ it[1].size() == 1 }
                     .map{ id, meta, path -> 
                        meta_notmerge=meta[0]
                        path_notmerge=path[0]
                        [meta_notmerge, path_notmerge]
                     }

    //ch_nocat.view()                 

    ch_sicat = INPUT_CHECK.out.reads
                     .map{ meta, path -> 
                        id=meta.subMap('id')
                        meta=meta
                        path=path
                        [id.id, meta, path]
                      }
                     .groupTuple()
                     .filter{ it[1].size() >= 2 } //filtra per numero di meta presenti dopo il tupla se hai due meta vuol dire che devi unire due campioni 
                     .map{ id, meta, path -> 
                        single = meta[0].subMap('single_end')
                        meta = meta[0]
                        //def flatPath = path.flatten()
                        [ meta , path[0], path[1] ]
                      }
                  

    CAT_FASTQ(
                        ch_sicat
                        ).reads.set { cat_fastq }


    //cat_fastq.view()
    //TODO authomatcly select the right ch_starter based on ss 
    //ch_starter = cat_fastq.mix(ch_nocat)
    //ch_starter = ch_nocat
    //ch_starter = cat_fastq

    //ch_nocat.view()

    ch_in = cat_fastq.mix(ch_nocat)
    //ch_in.view()

    ch_starter=ch_in
    
    //ch_starter = cat_fastq.ifEmpty(ch_in)
    
    //ch_starter = ch_nocat.ifEmpty(cat_fastq) 
       
    // if (Channel.exists('cat_fastq')) {

    //     ch_in = cat_fastq.mix(ch_nocat)
    //     ch_starter=ch_in
    //     ch_starter = cat_fastq.ifEmpty(ch_in) 

    // } else {
    //     ch_starter = ch_nocat
    // }

    //ch_starter.view()

    //
    // MODULE: Run FastQC
    //

    
    //ch_starter = merge_reads.mix(ch_nocat)

    FASTQC (
        //INPUT_CHECK.out.reads
        //cat_fastq
        ch_starter
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // ALLIGNMENT
    //
    // FASTQ_ALIGN_BWA(
    //     //cat_fastq,
    //     ch_starter,
    //     [[],params.index],
    //     false,
    //     [[],params.fasta]
    // )


    if (params.fillRE == false){
        BWA_MEM ( ch_starter, 
                    [[],params.index], 
                    [[],params.fasta], 
                    false 
                )
    }else{
       BWA_ALN (
          ch_starter,
          [[],params.fasta],
          //[,],
          [[],params.index]
       )

      //  BWA_SAMSE (

      //  )
      //TODO convert in bam
    }



    //
    //
    //

    //ch_bamin=FASTQ_ALIGN_BWA.out.bam
    ch_bamin=BWA_MEM.out.bam

    ch_chrom_sizes = CHROMSIZES ( params.fasta ).sizes

    PAIRTOOLS_PARSE(
        ch_bamin,
        params.fasta
    )

    PAIRTOOLS_SORT(
        PAIRTOOLS_PARSE.out.pairsam
    )

    PAIRTOOLS_DEDUP(
       PAIRTOOLS_SORT.out.sorted 
    )

    //ch_pairs=PAIRTOOLS_SORT.out.pairs

    PAIRTOOLS_SPLIT(
       PAIRTOOLS_DEDUP.out.pairs
    )

    PAIRTOOLS_QC (
        PAIRTOOLS_DEDUP.out.stat
    )

    BAM_SORT_STATS_SAMTOOLS(
        PAIRTOOLS_SPLIT.out.bam,
        [[],params.fasta]
    )

    //TODO make docker for dependencies
    //PAIRTOOLS_QC (
    //    PAIRTOOLS_DEDUP.out.stat
    //)

    //TODO  have commented this //args = task.attempt > 1 ? args.join(' -defects') : args  // Disable testing for defects 
    // find other solution

    if (params.exclude_lcextrap == false){
        PRESEQ_LCEXTRAP (
            BAM_SORT_STATS_SAMTOOLS.out.bam,
        )
    }

    //TODO add juicer
    // JUICERTOOLS(
    //     PAIRTOOLS_DEDUP.out.pairs,
    //     ch_chrom_sizes
    // )

    PAIRIX(
        PAIRTOOLS_DEDUP.out.pairs
    )

    PAIRIX.out.index.map{meta,path,index -> 
                    [meta,path,index,params.cool_bins]
                    }.set{ch_cloadin}
    
    COOLER_CLOAD(
        ch_cloadin,
        ch_chrom_sizes
    )


    if (params.skip_juicer == false){

        JUICER_PRE_50k(
            PAIRTOOLS_DEDUP.out.pairs,
            params.juicertool_location,
            ch_chrom_sizes,
            50000
        )

        JUICER_PRE_250k(
            PAIRTOOLS_DEDUP.out.pairs,
            params.juicertool_location,
            ch_chrom_sizes,
            250000
        )

        JUICER_PRE(
            PAIRTOOLS_DEDUP.out.pairs,
            params.juicertool_location,
            ch_chrom_sizes,
            params.cool_bins
        )

    }




    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMicroc.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMicroc.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
