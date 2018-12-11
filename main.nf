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
      --chr_file                    chromosome list to split bedfile

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
params.chr_file = false

fasta_file = check_ref_exist(params.fasta, 'fasta')
gtf_gff_file = check_ref_exist(params.gtf_gff, 'gtf/gff')
chr_file = check_ref_exist(params.chr_file, 'chr list')
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
if (!samtools_idx.exists()) {
    process samtools_index {
        tag "Build samtools fasta index on ${fasta_file.getName()}"

        publishDir outdir, mode: 'copy'

        input:
        file fasta from fasta_file

        output:
        file "${fasta}.fai" into genome_fai

        script:
        """
        samtools faidx ${fasta}
        """
    } 
} else {
    genome_fai = samtools_idx
}


//Prepare picard index
picard_idx = file("${fasta_file}.dict")
if (!picard_idx.exists()) {
    process picard_index {
        tag "Build picard fasta index on ${fasta_file.getName()}"

        publishDir outdir, mode: 'copy'

        input:
        file fasta from fasta_file

        output:
        file "${fasta.baseName}.dict" into genome_dict

        script:
        """
        java -jar /public/software/picard/2.18.17/picard.jar CreateSequenceDictionary \\
            R=${fasta} \\
            O="${fasta.baseName}.dict"
        """
    }
} else {
    genome_dict = picard_idx
}


// Prepare bwa index
bwa_idx_files = file("${outdir}/bwa/*")

if (bwa_idx_files.isEmpty()) {
    process bwa_index {
        tag "Build bwa index on ${fasta_file.getName()}"

        publishDir "${outdir}/bwa", mode: 'copy'

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
} else {
    bwa_idx = bwa_idx_files
}


// Prepare cds&exon bed
cds_exon_bed_files = file("${outdir}/*[cds,exon].merged.bed")

if (cds_exon_bed_files.isEmpty()) {
    process cds_exon_bed {
        tag "Exon & cds bed file from ${gtf.getName()}"

        publishDir outdir, mode: 'copy'

        when:
        cds_exon_bed_files.isEmpty()

        input:
        file gtf from gtf

        output:
        file "${gtf.baseName}.exon.merged.bed" into exon_bed
        file "${gtf.baseName}.cds.merged.bed" into cds_bed
        
        script:
        """
        python ${script_dir}/exon_cds_bed.py \\
            ${gtf.getName()}
        """
    }
} else {
    exon_bed = file("${outdir}/*.exon.merged.bed")
    cds_bed = file("${outdir}/*.cds.merged.bed")
}


// Split bed by chromosome list
process split_bed_by_chrom {

    publishDir "${outdir}", mode: 'copy'

    input:
    file bedfile from exon_bed
    file fai from genome_fai
    file chr_file from chr_file
    
    output:
    file "split_bed" into split_bed
    
    script:
    """
    python ${script_dir}/split_bed.py \\
        --input-file ${bedfile} \\
        --chr-list ${chr_file} \\
        --outdir ./split_bed/exon

    python ${script_dir}/split_bed.py \\
        --input-file ${fai} \\
        --chr-list ${chr_file} \\
        --fai True \\
        --outdir ./split_bed/genome
    """
}

// store index information int nextflow config
process make_config {
    tag "Generate config file: ${cfg_name}"

    publishDir outdir, mode: 'copy'

    input:
    file fasta from fasta_file
    file gtf from gtf
    file exon_bed from exon_bed
    file cds_bed from cds_bed
    file split_bed from split_bed
    file genome_dict from genome_dict
    file bwa_idx from bwa_idx

    output:
    file "${cfg_name}.nf.config"

    script:
    """
    echo "
    params {
        fasta = \'${fasta_path}/${fasta.getName()}\'
        gtf = \'${gtf_path}/${gtf.getName()}\'
        bwa_index = \'${outdir}/bwa/\'
        exon_bed = \'${outdir}/${gtf.baseName}.exon.bed\'
        cds_bed = \'${outdir}/${gtf.baseName}.cds.bed\'
        exon_split_bed = \'${outdir}/split_bed/exon\'
        genome_split_bed = \'${outdir}/split_bed/genome\'
    }
    " > ${cfg_name}.nf.config
    """
}