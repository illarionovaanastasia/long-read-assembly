#!/usr/bin/env nextflow
/*
========================================================================================
                         long-read-assembly
========================================================================================
 long-read-assembly Analysis Pipeline. Started 2020-07-15.
 #### Homepage / Documentation
 https://github.com/illarionovaanastasia/long-read-assembly
 #### Authors
 Anastasia Illarionova illarionovaanastasia <anastasia.illarionova@dzne.de> - https://github.com/illarionovaanastasia>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     long-read-assembly v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run illarionovaanastasia/long-read-assembly --Reads 'nano_reads.fastq.gz' --fasta reference.fasta --genomeSize 3300000000 -profile docker --sv_detection run 

    Required arguments
      --Reads                   Long reads
	  --genomeSize          Size of haploid genome in bp
      -profile                     Hardware config to use.

    References                 
      --fasta                     Path to Fasta reference

    Options:
      --lr_type                     Long read technology. One of 'nanopore' | 'pacbio' . Default: 'nanopore'
      --cont                        What container tech to use. Set to 'docker' if you use docker.
	  --sv_detection            Optional step to run SV detection on the de novo assembly if reference was provided. Default: 'false'

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
        .fromPath( params.Reads )
        .ifEmpty { exit 1, "Cannot find any long reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .into { reads_qc; reads_map_phasing; reads_map_polishing }


///*
// * Create a channel for reference assembly
// */
Channel
        .fromPath( params.fasta )
        .ifEmpty { exit 1, "Cannot find a reference genome assembly: ${params.fasta}\nNB: Path needs to be enclosed in quotes!" }
        .into { reference_map_phasing; reference_phasing; reference_scaffolding; sv_reference}


// Header log info
log.info "========================================="
log.info " long-read-assembly v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Long Reads']   = params.Reads
summary['Fasta Ref']    = params.fasta
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


 
/**
 * STEP 1 Reads QC
 * Only executed if docker is used
 */
process nanoqc {
    tag "${lreads.baseName}"
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




/**
 * STEP 2 FLYE Assembly
 */

if (params.genomeSize == 0){
  log.error "No genome size specified. Necessary for Flye assembly workflow"
  exit 1
}

// Map short reads to assembly with minimap2
process minimap_phasing {
  publishDir "${params.outdir}/minimap_phasing", mode: 'copy'
  
  input:
    file reference from reference_map_phasing
    file lreads from reads_map_phasing
  
  output:
    file "output.bam" into aligned_reads
	
  script:
  
  rtype = (params.lr_type == "nanopore") ? "map-ont" : "map-pb"
    """
	minimap2 -ax $rtype $reference $lreads |  samtools sort -o output.bam
	
  """
  
}

// Phase reads with MarginPhase
process marginphase {
  publishDir "${params.outdir}marginphase", mode: 'copy'
  
  input:
    file alignment from aligned_reads
    file reference from reference_phasing
  
  output:
    file "hap1.fastq" into hap1_reads
    file "hap2.fastq" into hap2_reads	
    file "*" into marginphase_results
  
  script:
  param = (params.lr_type == "nanopore") ? "params.nanopore.json" : "params.pacbio.json"
    """
  marginPhase $alignment $reference /marginPhase/params/$param
  samtools view -S -b output.1.sam > phasing.1.bam
  samtools view -S -b output.2.sam > phasing.2.bam
  samtools bam2fq phasing.1.bam > hap1.fastq
  samtools bam2fq phasing.2.bam > hap2.fastq

  
  """		
}


// Create assembly with FLYE
process flye {
  publishDir "${params.outdir}flye", mode: 'copy'
  
  input:
    file hap1 from hap1_reads
    file hap2 from hap2_reads	
  
  output:
    file "assembly1.fasta" into assembly_flye1, assembly_polishing_flye1
    file "assembly2.fasta" into assembly_flye2, assembly_polishing_flye2
  
  script:
    """
  flye --nano-raw $hap1 --genome-size $params.genomeSize --threads 60 --out-dir flye1
  mv flye1/assembly.fasta assembly1.fasta
     
  flye --nano-raw $hap2 --genome-size $params.genomeSize --threads 60 --out-dir flye2
  mv flye2/assembly.fasta assembly2.fasta
  
  """		
}


// Map long reads to assembly with minimap2
process minimap_polishing {
  publishDir "${params.outdir}/minimap_polishing", mode: 'copy'
  
  input:
    file lreads from reads_map_polishing
    file assembly1 from assembly_flye1
    file assembly2 from assembly_flye2
  
  output:
    file "output1.sorted.bam" into polishing_alignment1
    file "output2.sorted.bam" into polishing_alignment2
  
  script:
    rtype = (params.lr_type == "nanopore") ? "map-ont" : "map-pb"
    """
    minimap2 -ax $rtype $assembly1 $lreads |  samtools sort -o output1.sorted.bam
    minimap2 -ax $rtype $assembly2 $lreads |  samtools sort -o output2.sorted.bam
  """
  
}

  // Polish flye assembly with MarginPhase-HELEN pipeline
    process polishing {
        publishDir "${params.outdir}/polishing", mode: 'copy'

        input:
        file polishing_alignment1 from polishing_alignment1
        file polishing_alignment2 from polishing_alignment2
        file assembly_polishing_flye1 from assembly_polishing_flye1
        file assembly_polishing_flye2 from assembly_polishing_flye2

        output:
        file "helen_hap1.fasta" into ragoo_contigs1
        file "helen_hap2.fasta" into ragoo_contigs2

        script:
        """
		cp /helen/venv/bin/HELEN_r941_guppy344_human.pkl HELEN_r941_guppy344_human.pkl
		
		samtools index $polishing_alignment1
		mkdir marginpolish_images_hap1
        marginpolish $polishing_alignment1 $assembly_polishing_flye1 /helen/venv/bin/MP_r941_guppy344_human.json -t 60 -o marginpolish_images_hap1/ -f
        helen polish --image_dir marginpolish_images_hap1/ --model_path HELEN_r941_guppy344_human.pkl --batch_size 256 --num_workers 0 --threads 60 --output_dir helen_hap1/ --output_prefix polished_hap1
		mv helen_hap1/polished_hap1.fa helen_hap1.fasta
		
		samtools index $polishing_alignment2
		mkdir marginpolish_images_hap2
        marginpolish $polishing_alignment2 $assembly_polishing_flye2 /helen/venv/bin/MP_r941_guppy344_human.json -t 60 -o marginpolish_images_hap2/ -f
        helen polish --image_dir marginpolish_images_hap2/ --model_path HELEN_r941_guppy344_human.pkl --batch_size 256 --num_workers 0 --threads 60 --output_dir helen_hap2/ --output_prefix polished_hap2
		mv helen_hap2/polished_hap2.fa helen_hap2.fasta
        """
		

    }
	
	
   // Scaffold assembly with RaGOO
    process ragoo {
        publishDir "${params.outdir}/ragoo_scaffolds", mode: 'copy'

        input:
        file assembly1 from ragoo_contigs1
        file assembly2 from ragoo_contigs2
        file r_assembly from reference_scaffolding

        output:
        file "ragoo1.fasta" into ragoo_scaffolds_fasta1, mumandco_scaffolds_fasta1
        file "ragoo2.fasta" into ragoo_scaffolds_fasta2, mumandco_scaffolds_fasta2

        script:
        """
        ragoo.py $assembly1 $r_assembly 
		mv ragoo_output ragoo_output1
		cp ragoo_output1/ragoo.fasta ragoo1.fasta
		
        ragoo.py $assembly2 $r_assembly 
		mv ragoo_output ragoo_output2
		cp ragoo_output2/ragoo.fasta ragoo2.fasta
        """

    }
	

// Quast for flye pipeline
process quast_flye {
  publishDir "${params.outdir}/quast_results", mode: 'copy'
  
  input:
    file scaffolds1 from ragoo_scaffolds_fasta1
	file scaffolds2 from ragoo_scaffolds_fasta2
  
  output:
    file "*" into quast_results
  
  script:
    """
    quast  $scaffolds1
	quast  $scaffolds2

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
        file scaffolds1 from mumandco_scaffolds_fasta1
		file scaffolds2 from mumandco_scaffolds_fasta2
        file r_assembly from sv_reference

        output:
        file "*" into mumandco

        script:
		"""
        bash mumandco.sh -r $r_assembly -q $scaffolds1 -o hap1 -g $params.genomeSize
		bash mumandco.sh -r $r_assembly -q $scaffolds2 -o hap2 -g $params.genomeSize
		"""
	}
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {
    log.info "[long-read-assembly] Pipeline Complete"

}
