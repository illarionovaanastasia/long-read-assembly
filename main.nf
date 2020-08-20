#!/usr/bin/env nextflow
/*
========================================================================================
                         long-read-assembly
========================================================================================
 long-read-assembly Analysis Pipeline. Started 2020-07-15.
 #### Homepage / Documentation
 https://github.com/illarionovaanastasia/longread-haploid-assembly
 #### Authors
 Anastasia Illarionova illarionovaanastasia <avalillarionova@gmail.com> - https://github.com/illarionovaanastasia>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     long-read-assembly v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run illarionovaanastasia/longread-haploid-assembly --Reads 'nano_reads.{1,2}.fastq.gz' --fasta reference.fasta --genomeSize 3300000000 -profile docker --sv_detection run

    Required arguments
      --Reads                   Long reads from 1 and 2 haplotypes
      --genomeSize              Size of haploid genome in bp
      -profile                     Hardware config to use.

    References                 
      --fasta                     Path to Fasta reference

    Options:
      --lr_type                     Long read technology. One of 'nanopore' | 'pacbio' . Default: 'nanopore'
      --cont                        What container tech to use. Set to 'docker' if you use docker.
      --sv_detection                Optional step to run SV detection on the de novo assembly if reference was provided. Default: 'false'
      --qc                            optional step to run QC with NanoPlot

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                      Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --name                      Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.fasta = ""
params.sv_detection = false
params.qc = false
params.Reads = ""
params.email = false
params.plaintext_email = false
params.genomeSize = false
params.cont = 'docker'

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
    log.info "Reference genome was specified. It will be used for scaffolding."
	if(params.qc) log.info "QC step will be run."
	if(params.sv_detection) log.info "SV-detection step will be run."
} else {
    fasta = file("placeholder") // create placeholder
}


// Name specified by user
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


///*
// * Create a channel for input long read files
// */


Channel
    .fromPath( params.Reads)
    .ifEmpty { exit 1, "Cannot find any long reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
    .into { reads_qc; reads_flye; reads_map_polishing }


///*
// * Create a channel for reference assembly
// */
Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "Cannot find a reference genome assembly: ${params.fasta}\nNB: Path needs to be enclosed in quotes!" }
        .into { reference_scaffolding; sv_reference}







// Header log info
log.info "========================================="
log.info " long-read-assembly v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Long Reads']   = params.Reads
summary['Fasta Ref']    = params.fasta
summary['QC']    = params.qc
summary['SV-detection']    = params.sv_detection
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="






// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}



// WORKFLOW FOR UNPHASED READS


/**
 * STEP 1 Reads QC
 * Only executed if docker is used
 */
 if (params.qc == 'run') {
    process nanoqc {
        publishDir "${params.outdir}/nanoqc", mode: 'copy'

        input:
        file lreads from reads_qc

        output:
        file "*" into nanoqc_results

        script:
        ftype = (lreads.extension == "fasta" || lreads.extension == "fa") ? "--fasta" : "--fastq"
        """
        source activate nanoqc-env
        NanoPlot $ftype $lreads

        """
        }
}

/**
 * STEP 2 FLYE Assembly
 */

if (params.genomeSize == 0){
  log.error "No genome size specified. Necessary for Flye assembly workflow"
  exit 1
}

// Create assembly with FLYE
process flye {
  publishDir "${params.outdir}flye", mode: 'copy'
  
  input:
    file lreads from reads_flye	
  
  output:
    file "assembly.fasta" into assembly_flye, assembly_polishing_flye


  script:
    """
  flye --nano-raw $lreads --genome-size $params.genomeSize --threads 60 --out-dir flye
  mv flye/assembly.fasta assembly.fasta
  
  """		
}


// Map long reads to assembly with minimap2
process minimap_polishing {
  publishDir "${params.outdir}/minimap_polishing", mode: 'copy'
  
  input:
    file lreads from reads_map_polishing
    file assembly from assembly_flye
  
  output:
    file "output.sorted.bam" into polishing_alignment
  
  script:
    rtype = (params.lr_type == "nanopore") ? "map-ont" : "map-pb"
    """
    minimap2  -t 60 -ax $rtype $assembly $lreads |  samtools sort -@ 60 -o output.sorted.bam


  """
  
}

  // Polish flye assembly with MarginPhase-HELEN pipeline
    process polishing {
        publishDir "${params.outdir}/polishing", mode: 'copy'

        input:
        file polishing_alignment from polishing_alignment
        file assembly_polishing_flye from assembly_polishing_flye


        output:
        file "helen.fasta" into ragtag_contigs

        script:
        """
	cp /helen/venv/bin/HELEN_r941_guppy344_human.pkl HELEN_r941_guppy344_human.pkl
		
	samtools index $polishing_alignment
	mkdir marginpolish_images
        marginpolish $polishing_alignment $assembly_polishing_flye /helen/venv/bin/MP_r941_guppy344_human.json -t 60 -o marginpolish_images/ -f
        helen polish --image_dir marginpolish_images/ --model_path HELEN_r941_guppy344_human.pkl --batch_size 256 --num_workers 0 --threads 60 --output_dir helen/ --output_prefix polished
	mv helen/polished.fa helen.fasta 

        """
		

    }
	
	
   // Scaffold assembly with RagTag
    process ragtag {
        publishDir "${params.outdir}/ragtag_scaffolds", mode: 'copy'

        input:
        file assembly from ragtag_contigs
        file r_assembly from reference_scaffolding

        output:
        file "ragtag.fasta" into ragtag_scaffolds_fasta, mumandco_scaffolds_fasta


        script:
        """
        ragtag.py $assembly $r_assembly ragtag_output/query.corrected.fasta 
	cp ragtag_output/query.corrected.fasta ragtag.fasta
		
        """

    }
	

// Quast for flye pipeline
process quast_flye {
  publishDir "${params.outdir}/quast_results", mode: 'copy'
  
  input:
    file scaffolds from ragtag_scaffolds_fasta
  
  output:
    file "*" into quast_results
  
  script:
    """
    quast  $scaffolds -o quast

    """
}	


    

/*
 * Step 3 SV-detection
 */
if (params.sv_detection == 'run') {
    process mumandco {
        tag "$name"
        publishDir "${params.outdir}/mumandco", mode: 'copy'

        input:
        file scaffolds from mumandco_scaffolds_fasta
        file r_assembly from sv_reference

        output:
        file "*" into mumandco

        script:
		"""
        bash mumandco.sh -r $r_assembly -q $scaffolds -o ${name} -g $params.genomeSize

		"""
	}
}







/*
 * Completion e-mail notification
 */
workflow.onComplete {
    log.info "[long-read-assembly] Pipeline Complete"

}
