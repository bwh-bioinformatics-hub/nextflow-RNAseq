#!/usr/bin/env nextflow

def json_schema = "$workflow.projectDir/nextflow_schema.json"

if(params.phenotype){
    pheno_file = file(params.phenotype)
    ch_phenotype = examine_phenotype(pheno_file)
} else {
    ch_phenotype = Channel.empty()
}

if(params.trim_fastq){
    if(params.adapters){
    adapters = file(params.adapters, checkIfExists: true)
    if(!params.k && !params.ktrim || !params.k && params.ktrim || params.k && !params.ktrim){
        exit 1, "[nf-core/circrna] error: Adapter file provided for trimming but missing values for '--k' and/or '--ktrim'.Please provide values for '--k' and '--ktrim'.\n\nPlease check the parameter documentation online."
    }
    }
    if(params.trimq && !params.qtrim || !params.trimq && params.qtrim){
        exit 1, "[nf-core/circrna] error: Both '--trimq' and '--qtrim' are required to perform quality filtering - only one has been provided.\n\nPlease check the parameter documentation online."
    }
}

tool = params.tool ? params.tool.split(',').collect{it.trim().toLowerCase()} : []

// Get Input data - could be fastq or bam files
csv_file = file(params.input, checkIfExists: true)
ch_input = extract_data(csv_file)

(bam_input, fastq_input) = ch_input.into(2)

// Stage config files
ch_multiqc_config = file("$workflow.projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()

/*
================================================================================
                        PRINTING PARAMETER SUMMARY
================================================================================
*/

log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName

def summary = [:]
if (workflow.revision)        summary['Pipeline Release'] = workflow.revision
summary['Run Name']           = custom_runName ?: workflow.runName
if (workflow.containerEngine) summary['Container'] = "${workflow.containerEngine} - ${workflow.container}"
summary['Max Resources']      = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
summary['Config Files']       = workflow.configFiles.join(', ')
summary['Launch dir']         = workflow.launchDir
summary['Output dir']         = params.outdir
summary['Publish dir mode']   = params.publish_dir_mode
summary['Working dir']        = workflow.workDir
summary['Script dir']         = workflow.projectDir
summary['User']               = workflow.userName

summary['Input']              = params.input
summary['Input type']         = params.input_type
summary['BSJ filter']         = params.bsj_reads

summary['Genome']             = params.genome
if(params.fasta)              summary['Reference FASTA']   = params.fasta
if(params.gtf)                summary['Reference GTF']     = params.gtf
if(params.bowtie)             summary['Bowtie indices']    = params.bowtie
if(params.bowtie2)            summary['Bowtie2 indices']   = params.bowtie2
if(params.bwa)                summary['BWA indices']       = params.bwa
if(params.fasta_fai)          summary['SAMtools index']    = params.fasta_fai
if(params.hisat)              summary['HISAT2 indices']    = params.hisat2
if(params.star)               summary ['STAR indices']     = params.star
if(params.segemehl)           summary ['Segemehl index']   = params.segemehl

summary['Save reference files']              = params.save_reference
summary['Save QC intermediates']             = params.save_qc_intermediates
summary['Save RNASeq intermediates']         = params.save_rnaseq_intermediates
summary['Save Quantification intermediates'] = params.save_quantification_intermediates
summary['Save miRNA predictions']            = params.save_mirna_predictions

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-circrna-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/circrna Workflow Summary'
    section_href: 'https://github.com/nf-core/circrna'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
================================================================================
                            Stage Parameters
================================================================================
*/

params.fasta     = params.genome ? params.genomes[params.genome].fasta ?: false : false
params.gtf       = params.genome ? params.genomes[params.genome].gtf ?: false : false
params.bwa       = params.genome ? params.genomes[params.genome].bwa ?: false : false
params.star      = params.genome ? params.genomes[params.genome].star ?: false : false
params.species   = params.genome ? params.genomes[params.genome].species_id?: false : false

ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : 'null'
ch_gtf = params.gtf ? Channel.value(file(params.gtf)) : 'null'
ch_mature = params.mature ? Channel.value(file(params.mature)) : 'null'
ch_species = params.genome ? Channel.value(params.species) : 'null'

/*
================================================================================
                            BUILD INDICES
================================================================================
*/
process STAR_INDEX {
    tag "${fasta}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/${it}" : null }

    when:
    !params.star && params.fasta && params.gtf

    input:
    file(fasta) from ch_fasta
    file(gtf) from ch_gtf

    output:
    file("STARIndex") into star_built

    script:
    """
    mkdir -p STARIndex

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbOverhang ${params.sjdbOverhang} \\
        --sjdbGTFfile $gtf \\
        --genomeDir STARIndex/ \\
        --genomeFastaFiles $fasta
    """
}

ch_star = params.star ? Channel.value(file(params.star)) : star_built

/*
================================================================================
                            Misc circRNA Requirements
================================================================================
*/
process FILTER_GTF{
    tag"${gtf}"

    input:
    file(gtf) from ch_gtf

    output:
    file("filt.gtf") into ch_gtf_filtered

    script:
    """
    grep -vf ${workflow.projectDir}/conf/unwanted_biotypes.txt $gtf > filt.gtf
    """
}

process GENE_ANNOTATION{
    tag "${gtf}"
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_reference ? "reference_genome/${it}" : null }

    when:
    params.gtf

    input:
    file(gtf) from ch_gtf

    output:
    file("${gtf.baseName}.txt") into ch_gene_txt

    script:
    """
    gtfToGenePred -genePredExt -geneNameAsName2 ${gtf} ${gtf.baseName}.genepred
    perl -alne '\$"="\t";print "@F[11,0..9]"' ${gtf.baseName}.genepred > ${gtf.baseName}.txt
    """
}

ch_gene = params.circexplorer2_annotation ? Channel.value(file(params.circexplorer2_annotation)) : ch_gene_txt

/*
================================================================================
                            Stage Input Data
================================================================================
*/
process BAM_TO_FASTQ{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { params.save_qc_intermediates ? "quality_control/SamToFastq/${it}" : null }

    when:
    params.input_type == 'bam'

    input:
    tuple val(base), file(bam) from bam_input

    output:
    tuple val(base), file('*.fq.gz') into fastq_built

    script:
    """
    picard \\
        -Xmx${task.memory.toGiga()}g \\
        SamToFastq \\
        I=$bam \\
        F=${base}_R1.fq.gz \\
        F2=${base}_R2.fq.gz \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

if(params.input_type == 'bam'){
    (fastqc_reads, trimming_reads, raw_reads) = fastq_built.into(3)
}else if(params.input_type == 'fastq'){
    (fastqc_reads, trimming_reads, raw_reads) = fastq_input.into(3)
}

rmats_phenotype, dea_phenotype = ch_phenotype.into(2)

process FASTQC_RAW {
    tag "${base}"
    label 'process_medium'
    label 'py3'

    input:
    tuple val(base), file(fastq) from fastqc_reads

    output:
    file("*.{html,zip}") into fastqc_raw

    script:
    """
    fastqc -q $fastq --threads ${task.cpus}
    """
}


/*
================================================================================
                                    TRIM
================================================================================
*/
process TRIM_GALORE {
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "*.fq.gz",
        saveAs: { params.save_qc_intermediates ? "quality_control/trimgalore/${it}" : null }

    when:
    params.trim_fastq

    input:
    tuple val(base), file(fastq) from trimming_reads

    output:
    tuple val(base), file('*.trim.fq.gz') into trim_reads_ch, trim_reads_ch_fastqc
    //file(*) into trim_results

    script:
    // Check the configs in ./nextflow.config
    def c_r1   = params.clip_r1 > 0             ? "--clip_r1 ${params.clip_r1}"                         : ''
    def c_r2   = params.clip_r2 > 0             ? "--clip_r2 ${params.clip_r2}"                         : ''
    def tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    def tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

    // check if pair reads or not
    if(fastq[1]){
    """
    trim_galore \\
        --cores ${task.cpus} \\
        --gzip \\
        --fastqc \\
        --paired \\
        $c_r1 \\
        $c_r2 \\
        $tpc_r1 \\
        $tpc_r2 \\
        ${fastq[0]} \\
        ${fastq[1]}
    """
    }
    else{
    """
    trim_galore \\
        --cores ${task.cpus} \\
        --gzip \\
        --fastqc \\
        $c_r1 \\
        $tpc_r1 \\
        ${fastq[0]}
    """
    }
}

aligner_reads = params.trim_fastq ? trim_reads_ch : raw_reads

process FASTQC_TRIMED {
    tag "${base}"
    label 'process_low'

    when:
    params.trim_fastq

    input:
    tuple val(base), file(fastq) from trim_reads_ch_fastqc

    output:
    file ("*.{html,zip}") into fastqc_trimmed

    script:
    """
    fastqc -q $fastq --threads ${task.cpus}
    """
}

/*
================================================================================
                                STAR Alignment
================================================================================
*/
process STAR{
    tag "${base}"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "STAR/${it}" : null }

    input:
    tuple val(base), file(reads) from aligner_reads
    file(star_idx) from ch_star

    output:
    tuple val(base), file("${base}/${base}.Chimeric.out.junction") into circexplorer2_input
    tuple val(base), file("${base}/${base}.Aligned.sortedByCoord.out.bam") into star_bam_rmats, star_bam_stringtie
    tuple val(base), file("${base}") into circrna_finder_input

    script:
    def readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
    """
    mkdir -p ${base}

    STAR \\
        --twopassMode Basic \\
        --alignIntronMax ${params.alignIntronMax} \\
        --alignIntronMin ${params.alignIntronMin} \\
        --alignMatesGapMax ${params.alignMatesGapMax} \\
        --alignSJDBoverhangMin ${params.alignSJDBoverhangMin} \\
        --alignSJoverhangMin ${params.alignSJoverhangMin} \\
        --alignSoftClipAtReferenceEnds ${params.alignSoftClipAtReferenceEnds} \\
        --alignTranscriptsPerReadNmax ${params.alignTranscriptsPerReadNmax} \\
        --chimJunctionOverhangMin ${params.chimJunctionOverhangMin} \\
        --chimOutType Junctions SeparateSAMold \\
        --chimScoreMin ${params.chimScoreMin} \\
        --chimScoreSeparation ${params.chimScoreSeparation} \\
        --chimSegmentMin ${params.chimSegmentMin} \\
        --genomeDir ${star_idx} \\
        --limitSjdbInsertNsj ${params.limitSjdbInsertNsj} \\
        --outFileNamePrefix ${base}/${base}. \\
        --outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread} \\
        --outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax} \\
        --outFilterMultimapNmax ${params.outFilterMultimapNmax} \\
        --outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange} \\
        --outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread} \\
        --outFilterType BySJout \\
        --outReadsUnmapped None \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outBAMsortingBinsN 150 \\
        --outSJfilterOverhangMin ${params.outSJfilterOverhangMin} \\
        ${readFilesCommand} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --sjdbScore ${params.sjdbScore} \\
        --winAnchorMultimapNmax ${params.winAnchorMultimapNmax}
    """
}

/*
================================================================================
                                rMATs Analysis
================================================================================
*/
process RMATS_ANALYSIS{
    tag "${base}"
    label 'process_high'
    publishDir "${params.outdir}/rMATS", mode: params.publish_dir_mode, pattern: "rmats-results*"
    
    // when:
    // param.if_enable_rmats

    input:
    file(bam) from star_bam_rmats.collect()
    file(phenotype) from rmats_phenotype
    file(gtf) from ch_gtf

    script:
    """
    python ${workflow.projectDir}/scripts/prepare_rmats_input.py ${phenotype}
    python ${workflow.projectDir}/scripts/execute_rmats.py ${gtf} ${task.cpus}
    """
}

/*
================================================================================
                                CircRNA Analysis
================================================================================
*/
process CIRCEXPLORER2{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_quantification_intermediates ? "circrna_discovery/CIRCexplorer2/${it}" : null }

    input:
    tuple val(base), file(chimeric_reads) from circexplorer2_input
    file(fasta) from ch_fasta
    file(gene_annotation) from ch_gene

    output:
    tuple val(base), val("CIRCexplorer2"), file("${base}") into circexplorer2_intermediates
    tuple val(base), file("${base}_circexplorer2.bed") into circexplorer2_results
    tuple val(base), val("CIRCexplorer2"), file("${base}_circexplorer2_circs.bed") into circexplorer2_annotated

    script:
    """
    mkdir -p ${base}

    CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}/${base}.STAR.junction.bed

    CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}/${base}.STAR.junction.bed -o ${base}/${base}.txt

    awk '{if(\$13 >= ${params.bsj_reads}) print \$0}' ${base}/${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_circexplorer2.bed

    ## Re-work for Annotation
    awk -v OFS="\t" '{print \$1, \$2, \$3, \$1":"\$2"-"\$3":"\$4, \$5, \$4}' ${base}_circexplorer2.bed > ${base}_circexplorer2_circs.bed
    """
}

process COUNT_MATRIX_SINGLE{
    publishDir "${params.outdir}/circrna_discovery", pattern: "count_matrix.txt", mode: params.publish_dir_mode

    input:
    file(bed) from circexplorer2_results.collect()
    val(tool) from params.tool

    output:
    file("circRNA_matrix.txt") into circRNA_counts
    file("count_matrix.txt") into matrix

    script:
    """
    # Strip tool name from BED files (no consolidation prior to this step for 1 tool)
    for b in *.bed; do
        basename=\${b%".bed"};
        sample_name=\${basename%"_${tool}"};
        mv \$b \${sample_name}.bed
    done

    python ${workflow.projectDir}/scripts/circRNA_counts_matrix.py > circRNA_matrix.txt
    Rscript ${workflow.projectDir}/scripts/reformat_count_matrix.R
    """
}

/*
================================================================================
                            Differential Expression
================================================================================
*/
process STRINGTIE{
    tag "${base}"
    label 'process_medium'
    publishDir params.outdir, mode: params.publish_dir_mode, pattern: "${base}",
        saveAs: { params.save_rnaseq_intermediates ? "differential_expression/intermediates/StringTie/${it}" : null }

    input:
    tuple val(base), file(bam) from star_bam_stringtie
    file(gtf) from ch_gtf

    output:
    file("${base}") into stringtie_dir

    script:
    """
    mkdir ${base}/
    stringtie $bam -e -G $gtf -C ${base}/${base}_cov.gtf -p ${task.cpus} -o ${base}/${base}.gtf -A ${base}/${base}_genes.list
    """
}

gene_type_mapping = file("${workflow.projectDir}/assets/my_gene_type.txt")

process DEA{
    label 'process_medium'
    publishDir "${params.outdir}/differential_expression", pattern: "RNA-Seq", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "circRNA", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "boxplots", mode: params.publish_dir_mode
    publishDir "${params.outdir}/differential_expression", pattern: "combined_differential_expression.csv", mode: params.publish_dir_mode
    publishDir "${params.outdir}/quality_control", pattern: "DESeq2_QC", mode: params.publish_dir_mode

    input:
    file(gtf_dir) from stringtie_dir.collect()
    file(circ_matrix) from circRNA_counts
    file(phenotype) from dea_phenotype
    val(species) from ch_species
    file(circexplorer_results) from circexplorer2_intermediates.collect()

    output:
    file("RNA-Seq") into rnaseq_dir
    // file("circRNA") into circrna_dir

    script:
    """
    for i in \$(ls -d */); do sample=\${i%"/"}; file=\${sample}.gtf; touch samples.txt; printf "\$sample\t\${i}\${file}\n" >> samples.txt; done

    prepDE.py -i samples.txt

    ## prepDE && circRNA counts headers are sorted where uppercase preceedes lowercase i.e Z before a
    ## reformat the phenotype file to match the order of the samples.
    head -n 1 $phenotype > header
    tail -n +2 $phenotype | LC_COLLATE=C sort > sorted_pheno
    cat header sorted_pheno > tmp && rm phenotype.csv && mv tmp phenotype.csv

    Rscript ${workflow.projectDir}/scripts/DEA.R gene_count_matrix.csv $phenotype $circ_matrix $species ${workflow.projectDir}/conf/ensemblDatabase_map.txt

    mv gene_count_matrix.csv RNA-Seq
    mv transcript_count_matrix.csv RNA-Seq

    ## webcrawler is disabaled in default
    python ${workflow.projectDir}/scripts/combination.py ${gene_type_mapping} disable_webcrawler
    """
}

/*
================================================================================
                            Pathway Enrichment
================================================================================
*/
pathway_enrich_script = file("${workflow.projectDir}/scripts/pathway_enrich.R")
process PATHWAY_ENRICHMENT{
    tag "${base}"
    label 'process_high'
    publishDir "${params.outdir}/pathway_enrichment", pattern: "results", mode: params.publish_dir_mode

    input:
    file(rnaseq) from rnaseq_dir
    

    script:
    """
    python ${workflow.projectDir}/scripts/trans_DE2csv.py
    python ${workflow.projectDir}/scripts/execute_pathway_enrichment.py
    """
}

/*
================================================================================
                                    MultiQC
================================================================================
*/
process MULTIQC{
    label 'process_low'
    publishDir "${params.outdir}/quality_control/MultiQC", mode: params.publish_dir_mode,
        pattern: "*.html"

    input:
    file(raw_fastqc) from fastqc_raw.collect().ifEmpty([])
    file(trim_fastqc) from fastqc_trimmed.collect().ifEmpty([])
    // file(trimgalore_stats) from trim_results.collect().ifEmpty([])
    file(multiqc_config) from ch_multiqc_config
    file(mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    // file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file("*.html") into multiqc_out

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc . -f $rtitle $rfilename $custom_config_file
    """
}

/*
================================================================================
                            Auxiliary functions
================================================================================
*/

// Check integer
def isValidInteger(value){
    value instanceof Integer
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}


// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "[nf-core] error:  Invalid CSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "[nf-core] error: Cannot find supplied FASTQ or BAM input file. If using input method CSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}"
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Read input files from input CSV
def extract_data(csvFile){
    Channel
        .fromPath(csvFile)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_keys = ["Sample_ID", "Read1", "Read2", "Bam"]
        if(!row.keySet().containsAll(expected_keys)) exit 1, "[nf-core] error: Invalid CSV input - malformed column names. Please use the column names 'Sample_ID', 'Read1', 'Read2', 'Bam'."

        checkNumberOfItem(row, 4)

        def samples = row.Sample_ID
        def read1 = row.Read1.matches('NA') ? 'NA' : return_file(row.Read1)
        def read2 = row.Read2.matches('NA') ? 'NA' : return_file(row.Read2)
        def bam = row.Bam.matches('NA') ? 'NA' : return_file(row.Bam)

        if(samples == '' || read1 == '' || read2 == '' || bam == '') exit 1, "[nf-core] error: a field does not contain any information. Please check your CSV file"
        if(read1.matches('NA') && read2.matches('NA') && bam.matches('NA')) exit 1, "[nf-core] error: A row in your CSV file appears to have missing information."
        if( !read1.matches('NA') && !has_extension(read1, "fastq.gz") && !has_extension(read1, "fq.gz") && !has_extension(read1, "fastq") && !has_extension(read1, "fq")) exit 1, "[nf-core/circrna] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
        if( !read2.matches('NA') && !has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz") && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "[nf-core/circrna] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
        if( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core] error: A specified BAM file either has a non-recognizable extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"

        // output tuple mimicking fromFilePairs if FASTQ provided, else tuple for BAM
        if(bam.matches('NA')){
            if(read2.matches('NA')){
                [ samples, read1 ]
            }else{
                [ samples, [read1, read2] ]
            }
        }else{
            [ samples, bam ]
        }

        }
}

// Check input phenotype file
def examine_phenotype(pheno){

  Channel
        .fromPath(pheno)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_cols = ['condition']

        if (!row.keySet().containsAll(expected_cols)) exit 1, "[nf-core/circrna] error: 'condition' is not a column name in the phenotype file.\n\nThe response variable must be named 'condition', please refer to the usage documentation online"

        def condition  = row.condition.matches('NA') ? 'NA' : row.condition

        if(condition == '') exit 1, "[nf-core/circrna] error: Invalid phenotype file, condition column contains empty cells."
        if(condition.matches('NA')) exit 1, "[nf-core/circrna] error: NA value in phenotype condition column."

        }
        .toList()

        return Channel.value(file(pheno))
}

/*
================================================================================
                            nf-core functions
================================================================================
*/

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/circrna v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

// handle multiqc_channels
//if(params.trim_fastq == false){
//    ch_multiqc_report = multiqc_trim_out
//}else{
//    ch_multiqc_report = multiqc_raw_out
//}

// Completion e-mail notification
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/circrna] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/circrna] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/circrna] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/circrna] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$workflow.projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$workflow.projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$workflow.projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$workflow.projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/circrna] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
            mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/circrna] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset  = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/circrna]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/circrna]${c_red} Pipeline completed with errors${c_reset}-"
    }
}
