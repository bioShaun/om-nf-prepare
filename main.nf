#!/usr/bin/env nextflow

def helpMessage() {
    log.info """

    A Pipeline to Prepare Reseq analysis refer

    Usage:

    ==================================================================

    nextflow run omReseq-prepare-nf --fasta genome.fa --gtf genome.gtf

    ==================================================================

    References
      --fasta                       Path to Fasta reference
      --gtf_gff                     Path to reference GTF file 

    Other options:
      --cfg_name                    Name of output configure file
      --outdir                      Path to store index files, default is directory store fasta

    """.stripIndent()
}

// workflow internal path&files
script_dir = file("$baseDir/script/")

/*
 * SET UP CONFIGURATION VARIABLES
 */


def check_ref_exist = {file_path, file_type ->
    if (file_path) {
        file_path = file(file_path)
        if( !file_path.exists() ) exit 1, "${file_type} file not found: ${file_path}"
        return file_path
    } else {
        exit 1, "No reference genome ${file_type} specified!"
    }
}


// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// check fasta/gtf file existence
params.fasta = false
params.gtf_gff = false
params.outdir = false
params.cfg_name = false

fasta_file = check_ref_exist(params.fasta, 'fasta')
gtf_gff_file = check_ref_exist(params.gtf_gff, 'gtf/gff')
fasta_path = fasta_file.getParent() 

//check output path
if (!params.outdir) {
    outdir = fasta_path
} else {
    outdir = params.outdir
}

// output config name
if (params.cfg_name) {
    cfg_name = params.cfg_name
} else {
    cfg_name = fasta_file.baseName
}

// gtf to gff
if (gtf_gff_file.getExtension() != 'gtf') {
    gtf_path = outdir
    process convertGFFtoGTF {
        tag "Convert gff3 [${gtf_gff_file.getName()}] to gtf"

        publishDir outdir, mode: 'copy'

        input:
        file gtf_gff from gtf_gff_file

        output:
        file "${gtf_gff_file.baseName}.gtf" into gtf

        script:
        """
        gffread $gtf_gff_file -T -o ${gtf_gff_file.baseName}.gtf
        """
    }
} else {
    gtf_path = gtf_gff_file.getParent()
    gtf = gtf_gff_file
}


//Prepare samtools index
samtools_idx = file("${fasta_file}.fai")

process samtools_index {
    tag "Build samtools fasta index on ${fasta_file.getName()}"

    publishDir outdir, mode: 'copy'

    when:
    !samtools_idx.exists()

    input:
    file fasta from fasta_file

    output:
    file "${fasta}.fai"

    script:
    """
    samtools faidx ${fasta}
    """
}

//Prepare picard index
picard_idx = file("${fasta_file}.dict")

process picard_index {
    tag "Build picard fasta index on ${fasta_file.getName()}"

    publishDir outdir, mode: 'copy'

    when:
    !picard_idx.exists()

    input:
    file fasta from fasta_file

    output:
    file "${fasta.baseName}.dict"

    script:
    """
    java -jar /public/software/picard/2.18.17/picard.jar CreateSequenceDictionary \
        R=${fasta} \
        O="${fasta.baseName}.dict"
    """
}

// Prepare bwa index
bwa_idx_files = file("${outdir}/bwa/*")

process bwa_index {
    tag "Build bwa index on ${fasta_file.getName()}"

    publishDir "${outdir}/bwa", mode: 'copy'

    when:
    bwa_idx_files.isEmpty()

    input:
    file fasta from fasta_file

    output:
    file "${fasta.getName()}*" into bwa_idx

    script:    
    """
    bwa index ${fasta}
    mkdir -p ${outdir}/bwa/
    ln -s ${fasta_path}/${fasta.getName()} ${outdir}/bwa/${fasta.getName()}
    """
}

// Prepare cds&exon bed
cds_exon_bed_files = file("${outdir}/*[cds,exon].bed")

process cds_exon_bed {
    tag "Exon & cds bed file from ${gtf.getName()}"

    publishDir outdir, mode: 'copy'

    when:
    cds_exon_bed_files.isEmpty()

    input:
    file gtf from gtf

    output:
    file "*.merged.bed" into cds_exon_bed
    
    script:
    """
    python ${script_dir}/exon_cds_bed.py \
        ${gtf.getName()}
    """
}

// store index information int nextflow config
process make_config {
    tag "Generate config file: ${cfg_name}"

    publishDir outdir, mode: 'copy'

    input:
    file fasta from fasta_file
    file gtf from gtf

    output:
    file "${cfg_name}.nf.config"

    script:
    """
    echo "
    params {
        fasta = ${fasta_path}/${fasta.getName()}
        gtf = ${gtf_path}/${gtf.getName()}
        bwa_index = "${outdir}/bwa/${fasta.getName()}"
        exon_bed = ${outdir}/${gtf.baseName}.exon.bed
        cds_bed = ${outdir}/${gtf.baseName}.cds.bed
    }
    " > ${cfg_name}.nf.config
    """
}