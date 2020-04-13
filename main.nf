#!/usr/bin/env nextflow
nextflow.preview.dsl=2
initParamsToDefaults(params)
if (params.help) exit 0, helpMessage()
_PLATFORM = "ILLUMINA"
_THREADS  = 32
templateDir = "${workflow.projectDir}/lib/bash_templates"
// processParams()
/* Process the parameters and set the environemnt */
params.name = 'Layer Lab DNA Seq Analysis Pipeline'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>
params.contactMail ?: 'mahmood.javaid@colorado.edu'
if (!params.sampleTsv && !params.reads){
    exit 1, 'You must sepcify either a sample tsv file or the reads, see --help'
}

if (!params.conditionsTsv){
    exit 1, 'You must sepcify a tsv showing experimental conditions see --help'
}
// Check parameters for profile fiji
// If they are not provided at the command line, set them for fiji from the env section of the config
if (workflow.profile == 'fiji'){
  resultsDir = params.resultsDir ?: "$fiji_results_dir"
  dbsnp = params.dbsnp ?: "$fiji_dbsnp"
  ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$fiji_ensembl_gene_annotation"
  known_indels = params.knownIndels ?: "$fiji_known_indels"
  ref_fasta = params.refFasta ?: "$fiji_ref_fasta"
  accessibleRegions = params.accessibleRegions?: "$fiji_accessible_regions"
  ref_dict = params.refDict ?: "$fiji_ref_dict"
  afOnlyGnomad = params.afOnlyGnomad ?: "$fiji_af_only_gnomad"
  afOnlySNPOnlyGnomad = params.afOnlySNPOnlyGnomad ?: "$fiji_af_only_snp_only_gnomad"
  // gnomad_genome = params.gnomadGenome?: "$fiji_gnomad_genome"
  gnomadRef = params.gnomadRef ?: "$fiji_gnomad_ref"
  homo_sapiens_genes = params.homoSapiensGenes ?: "$fiji_homo_sapiens_genes"
  homo_sapiens_exons = params.homoSapiensExones ?: "$fiji_homo_sapiens_exons"
  cosmicCodingVCF = params.cosmicCodingVcf ?: "$fiji_cosmic_coding"
  cosmicGenes = params.cosmicGenes ?: "$fiji_cosmic_genes"
  intervals_list = params.intervalsList ?: "$fiji_intervals_list"
  intervalsListUnAnn = params.intervalsList ?: "$fiji_intervals_list_unannotated"
  cadd = params.cadd ?: "$fiji_cadd"
  mpc= params.mpc ?: "$fiji_mpc"

  
}
// Check parameters for mendel
// if (workflow.profile == 'mendel'){
//   resultsDir = params.resultsDir ?: "$mendel_results_dir"
//   dbsnp = params.dbsnp ?: "$mendel_dbsnp"
//   ensemblGeneAnnotation = params.ensemblGeneAnnotation ?: "$mendel_ensembl_gene_annotation"
//   knownIndels = params.knownIndels ?: "$mendel_known_indels"
//   refFasta = params.refFasta ?: "$mendel_ref_fasta"
// }
// Finally check if all the required parameters have either been provided via the commandline or
// by specifying a particular profile such as fiji or mendel
if (!resultsDir) exit 1, 'Specify the directory to store analysis results, see --help'  
if (!dbsnp) exit 1, 'Specify the dbsnp (*.vcf.gz) file, see --help'  
if (!ensemblGeneAnnotation) exit 1, 'Specify the ensemblGeneAnnotation (*.gff3.gz) file, see --help'  
if (!known_indels) exit 1, 'Specify the known indels (*.vcf.gz) file, see --help'  
if (!ref_fasta) exit 1, 'Specify the refenrece fasta (*.fasta.gz) file, see --help'  
if (!afOnlyGnomad) exit 1, 'Specify the af_only_gnomad (*.vcf.gz) file, see --help'  
// if (!intervals_list) exit 1, 'Specify the intervals list (*.vcf.gz) file, see --help'  

if (!params.runName) exit 1, 'Specify a run name, see --help'
runNameWithoutSpaces = params.runName.replaceAll('\\W','_')
OUT_DIR = "${resultsDir}/${runNameWithoutSpaces}"
log.info "Output dir: ${OUT_DIR}"

if (!workflow.commandLine.contains('-resume') && file(OUT_DIR).exists())
  exit 1, "${OUT_DIR} already exists, specify another results dir or run name"

/* End of parameter handling section */

// Handle input files
inputFiles = Channel.empty()
tsvPath = ''
if (params.sampleTsv) tsvPath = params.sampleTsv

if (tsvPath){
    tsvFile = file(tsvPath)
    if (!tsvFile.exists()) exit 1, "${tsvPath} does not exists!"
    inputFiles = LLabUtils.extractSamples(tsvFile)
}else{
  // Define channel for reading file pairs
  Channel
      .fromFilePairs( params.reads, size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}"}
      .set {inputFiles}
}
list_of_controls = LLabUtils.extractControls(file(params.conditionsTsv))
list_of_cases = LLabUtils.extractCases(file(params.conditionsTsv))

if (params.withoutGatk & params.withoutSavvycnv & params.withoutCnvkit){
  println ("You have disablled all three callers (gatk, savvy, cnvkit), only the allignment will run")
}
else{
  println ("Following copy number callers will run: ")
  if (!params.withoutGatk){
    println ("GATK (mutect2 based)")
  }

  if (!params.withoutSavvycnv){
    println ("SavvyCNV")
  }

  if (!params.withoutCnvkit){
    println ("CNVKit")
  }
}

// sample_name = 'KU1919P'
// if (list_of_normals.find{it == sample_name}){
//   println("Found $sample_name")
// }else{
//   println("Did not find $sample_name")
// }

// file(params.conditionsTsv)
//     .readLines()
//     .each { println it }
// ch_conditions = LLabUtils.extractConditions(file(params.conditionsTsv))

// ch_conditions = LLabUtils.extractConditions(file(params.conditionsTsv))
// CRAMS_DIR = "${OUT_DIR}/aligned/crams"
// BAMS_DIR = "${OUT_DIR}/aligned/bams"
// RECAL_TABLES_DIR="${OUT_DIR}/misc/recal_table"
// ch_crams = Channel.fromPath( "${CRAMS_DIR}/*.cram" )
// ch_crams_index = Channel.fromPath( "${CRAMS_DIR}/*.cram.crai" )
// ch_bams = Channel.fromFilePairs( "${BAMS_DIR}/*.{bam,bam.bai}" )
// ch_crams.
// view{println it}
// ch_crams = 
//   Channel
//     .fromFilePairs("${CRAMS_DIR}/*.{cram,cram.crai}", size: -1) 
//     { file -> file.simpleName }
//     // .println { ext, files ->  "Files with the extension $ext are $files" }

// ch_bams = 
//   Channel
//     .fromFilePairs("${BAMS_DIR}/*.{bam,bam.bai}", size: -1) 
//     { file -> file.simpleName }
// PON_PATH='/Shares/layer_shared/projects/sequence_analysis/cancer_center_hg38/misc/pon/pon.vcf.gz'
// ch_pon = Channel.fromPath(PON_PATH)
// Generate a channel holding chromosome list
listOfChromosoms = []
if (workflow.profile == 'fiji_hg37'){
    listOfChromosoms = LLabUtils.getChrmList()
}else{
   listOfChromosoms = LLabUtils.getChrmListHg38()
}

Channel.from(LLabUtils.getChrmList())
    .set{chrmListCh}

// CRAMS_DIR_ABS_PATH = "${outDir}/align"

workflow{
    MapReads(inputFiles)
    CramToBam(MapReads.out[0])
    if (params.withoutGatk == false)
    {
      // GATK portion
      PreprocessIntervals()      
      CollectReadCounts(MapReads.out,
                        PreprocessIntervals.out)
      CreateReadCountSomaticPON(
        CollectReadCounts.out.collect()
      )
      
      DenoiseReadCounts(
        CollectReadCounts.out,
        CreateReadCountSomaticPON.out
      )
      PlotDenoisedCopyRatios(DenoiseReadCounts.out[0],
                            DenoiseReadCounts.out[1])

    

      RunModelSegments(DenoiseReadCounts.out[1])
      PlotModeledSegments(DenoiseReadCounts.out[1].collect(),
                      RunModelSegments.out[2])
      CallCopyRatioSegments(RunModelSegments.out[1])
      GetGeneSummaryGATK(CallCopyRatioSegments.out)
    }
    if (params.withoutSavvycnv == false)
    {
      // SavvyCNV related Summaries
      GenSavvyCNVCoverageSummary(CramToBam.out)
      RunSavvyCNV(GenSavvyCNVCoverageSummary.out.collect())
    }
      // TestProc(
      //   CramToBam.out[0].collect(),
      //   CramToBam.out[1].collect()
      // )
        
    if (params.withoutCnvkit == false)
    {
      // CNVKIT related portion
      RunCNVKitBatchMode(
        CramToBam.out[0].collect(),
        CramToBam.out[1].collect()
      )

    }
} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

def initParamsToDefaults(params){
  def null_path = ''
  params.reads = null_path
  params.runName = ''
  params.contactMail=''
  // params.resultsDir = null_path
  params.dbsnp = null_path
  params.ensemblGeneAnnotation = null_path
  params.refFasta = null_path
  params.knownIndels = null_path

  params.help = false  
  params.skipQc = false
  params.skipFastQc = false
  params.withoutGatk = false
  params.withoutSavvycnv = false
  params.withoutCnvkit = false
}

/* Helper functions */

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}
def layerLabMessage() {
  log.info "Workflow For Somatic Variations ~ ${workflow.manifest.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def helpMessage() {
  // Display help message
  this.layerLabMessage()
  runStrGitHub = "nextflow run javaidm/layer_lab_somatic_sv "
  runStr = "nextflow run main.nf "
  log.info "    Usage:"
  log.info "       $runStr --params-file <params.yaml>"
  log.info " * All the params below could be set (in Camel Case) in a yaml file (see the params.yaml as an exampl3) *"
  log.info "    --sample-tsv <file.tsv>"
  log.info "       Specify a TSV file containing sample_id, path_to_read1, path_to_read2."
  log.info "    --results-dir <directory>"
  log.info "       Specify a directory to hold the analysis results."
  log.info "    --run-name"
  log.info "       Specify a run name, results will be stored under results-dir/run-name"
  log.info "    --conditions-tsv <conditions.tsv>"
  log.info "       Specifying conditions of the experiment"
  log.info "    --without-gatk "
  log.info "       Run the pipeline without the GATK cnv workflow (mutect2 based)"
  log.info "    --without-savvycnv "
  log.info "       Run the pipeline without the SavvyCNV cnv caller"
  log.info "    --without-cnvkit "
  log.info "       Run the pipeline without the CNVKit cnv caller"
  log.info "    --help"
  log.info "       you're reading it."
  
}

/* Processes */

process MapReads {
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/align/crams/" , mode: 'copy', overwrite: false
    input:
    tuple val(sample), file(reads)

    output:
    file(out_file)
    file("${out_file}.crai")

    script:
    r1 = reads[0]
    r2 = reads[1]
    
    script:
    out_file = "${sample}.cram"
   """
   bwa mem \
    -t $_THREADS -R "@RG\tID:$sample\tSM:$sample\tPL:$_PLATFORM\tPU:$sample\tLB:$sample" $ref_fasta $r1 $r2 \
    | samblaster \
    | samtools sort --output-fmt-option seqs_per_slice=4000 -O CRAM --reference $ref_fasta -m 18G -@ 6 /dev/stdin -o $out_file \
    && samtools index $out_file
    
   """
}

process CramToBam{
    echo true
    tag "$sample"
    publishDir "${OUT_DIR}/align/bams" , mode: 'copy', overwrite: false
    input:
    file(cram)

    output:
    file(out_file)
    file("${out_file}.bai")

    script:
    sample = "${cram.simpleName}"
    out_file = "${sample}.bam"
    
    script:
   """
   samtools view -T $ref_fasta -b -o $out_file $cram \
   && samtools index $out_file
   """
}


process PreprocessIntervals {
    publishDir "${OUT_DIR}/misc/preprocessed_intervals", mode: 'copy'
    // cache false
    input:
    // file (interval_list)
    
    output:
    file(out_file)
    script:
    out_file = 'preprocessed.interval_list'
    """
    gatk PreprocessIntervals \
        -L $intervals_list \
        -R $ref_fasta \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $out_file
    """
}

process CollectReadCounts {
    echo true
    tag "$sample"
    // cache false
    publishDir "${OUT_DIR}/misc/read_counts", mode: 'copy'

    input:
    file(cram)
    file(cram_index)
    file(preprocessed_intervals)

    output:
    file(out_file)

    script:
    // cram = cram_index_pair[0]
    // cram_index = cram_index_pair[1]
    sample = cram.simpleName
    out_file = "${sample}.counts.hdf5"
    
    
    """
    gatk CollectReadCounts \
        -I $cram \
        -R $ref_fasta \
        -L $preprocessed_intervals \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $out_file
    """
}

process CreateReadCountSomaticPON {
    echo true
    // tag "$sample"
    
    publishDir "${OUT_DIR}/misc/read_count_somatic_pon", mode: 'copy'
    
    input:
    // file (read_count_hdf5s:'all_read_counts/*')
    file(read_count_hdf5s)
    
    // file(this_read)

    output:
    file(out_file)

    script:
    // sample = this_read.simpleName
    out_file = "cnv.pon.hdf5"
    params_str = ''

    // Only get the normal samples
    read_count_hdf5s.each{
        sample = it.simpleName
        // check if this is one of the controls/parental samples
        if (LLabUtils.sampleInList(sample, list_of_controls)){
          params_str = "${params_str} -I ${it}"
        }
    }

    
    """
    gatk CreateReadCountPanelOfNormals \
        $params_str \
        -O $out_file
    """
}

process DenoiseReadCounts {
    echo true
    tag "$sample"
    
    publishDir "${OUT_DIR}/misc/denoised_read_counts", mode: 'copy'
    
    input:
    file(sample_counts_hdf5)
    file(read_count_somatic_pon)

    output:
    file(std_copy_ratio)
    file(denoised_copy_ratio)

    script:
    sample = sample_counts_hdf5.simpleName
    std_copy_ratio = "${sample}.standardizedCR.tsv"
    denoised_copy_ratio = "${sample}.denoisedCR.tsv"
    
    """
    
    gatk DenoiseReadCounts \
        -I $sample_counts_hdf5 \
        --count-panel-of-normals $read_count_somatic_pon \
        --standardized-copy-ratios $std_copy_ratio \
        --denoised-copy-ratios $denoised_copy_ratio
    """
}

process PlotDenoisedCopyRatios {
    echo true
    tag "$sample"
    
    publishDir "${OUT_DIR}/misc/denoisedCR_plots", mode: 'copy'
    
    input:
    file(std_copy_ratio)
    file(denoised_copy_ratio)

    output:
    file(out_dir)

    script:
    sample = std_copy_ratio.simpleName
    out_dir = "${sample}"
    
    """
    mkdir $out_dir
    gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios $std_copy_ratio \
        --denoised-copy-ratios $denoised_copy_ratio \
        --sequence-dictionary $ref_dict \
        --output-prefix $sample \
        -O $out_dir
    """
}
process RunModelSegments {
    echo true
    tag "$sample"
    
    publishDir "${OUT_DIR}/misc/modeled_segments", mode: 'copy', overwrite: true
    
    input:
    file(case_denoisedCR)

    output:
    file(out_dir)
    file("${out_dir}/${sample}.cr.seg")
    file("${out_dir}/${sample}.modelFinal.seg")

    script:
    sample = case_denoisedCR.simpleName
    out_dir = "${sample}_modeled_segments"
    """
    mkdir $out_dir
    gatk ModelSegments \
        --denoised-copy-ratios $case_denoisedCR \
        --output-prefix $sample \
        -O $out_dir
    """
}
process PlotModeledSegments {
    echo true
    tag "$sample"
    
    publishDir "${OUT_DIR}/misc/modeled_segments_plots", mode: 'copy', overwrite: true
    
    input:
    file("denoised_crs/*")
    file(model_final_seg)

    output:
    file(out_dir)

    script:
    sample = model_final_seg.simpleName
    sample_denoised_cr = "denoised_crs/${sample}.denoisedCR.tsv"
    out_dir = "${sample}"
    
    """
    mkdir $out_dir
    gatk PlotModeledSegments \
        --denoised-copy-ratios $sample_denoised_cr \
        --segments $model_final_seg \
        --sequence-dictionary $ref_dict \
        --output-prefix $sample \
        -O $out_dir
    """
}

process CallCopyRatioSegments {
    echo true
    tag "$sample"
    
    publishDir "$publish_dir", mode: 'copy', overwrite: true
    
    input:
    file(cr_seg)
    
    output:
    file(out_file)
    // file(publish_dir)

    script:
    sample = cr_seg.simpleName
    publish_dir = "${OUT_DIR}/misc/called_cr_segments"
    out_file = "${sample}.called.seg"
    
    """
    gatk CallCopyRatioSegments \
        -I $cr_seg \
        -O $out_file
    """
}
process GetGeneSummaryGATK{
    
    echo true
    publishDir "${OUT_DIR}/misc/gene_summary_gatk", mode: 'copy', overwrite: true

    input:
    file(called_seg)
    // file(bcf_index)

    output:
    file "$out_file"

    script:
    out_file = "${called_seg.simpleName}_gatk_gene_summary.txt"    
    
    """
    cat $called_seg | grep '^chr[1-9XY]' > ${called_seg}.modified
    bgzip ${called_seg}.modified
    tabix -s 1 -b 2 -e 3 ${called_seg}.modified.gz
    export f=${called_seg}.modified.gz  
    cat $homo_sapiens_genes \
    | gargs -p 8 'echo -en "{3}\t"; tabix \$f {0}:{1}-{2}' \
    |awk 'BEGIN {OFS = FS} {if (\$NF != 0) print \$1,\$NF}' > ${out_file}
    """
}
/* SavvyCNV Related processes */

process GenSavvyCNVCoverageSummary {
    echo true
    tag "$sample"
    // cache false
    publishDir "$publish_dir", mode: 'copy', overwrite: true
    
    input:
    file(bam)
    file(bam_index)
    
    output:
    file(out_file)
    // file(publish_dir)

    script:
    sample = bam.simpleName
    publish_dir = "${OUT_DIR}/misc/SavvyCNV_coverage_summaries"
    out_file = "${sample}.coverageBinner"
    
    """
    java -Xmx1g CoverageBinner $bam >$out_file 
    """
}

process RunSavvyCNV {
    echo true
    // cache false
    // tag "$sample"
    
    publishDir "$publish_dir", mode: 'copy', overwrite: true
    
    input:
    file("*")
    
    output:
    file("cnv_list.csv")
    file("log_messages.txt")
    file("*.pdf")
    // file()

    script:
    // sample = bam.simpleName
    publish_dir = "${OUT_DIR}/misc/SavvyCNV_calls"
    // out_file = "${sample}.coverageBinner"
    chunk_size = 200000
    
    """
    java -Xmx30g SavvyCNV -a -d $chunk_size *.coverageBinner >cnv_list.csv 2>log_messages.txt
    """
}
// CNVKIT related
process RunCNVKitBatchMode{
  echo true
  // cache false
  publishDir "$publish_dir", mode: 'copy', overwrite: true
  
  input:
    file("*")
    file("*")
  
  output:
  file('results')

  script:
  // Create a string from the list of controls and cases
  case_bams='' 
  control_bams = ''
  list_of_cases.each{ case_bams += "${it}.bam "}
  list_of_controls.each{ control_bams += "${it}.bam "}

  publish_dir = "${OUT_DIR}/misc/CNVKit_calls"
  """
  ls -al .
  cnvkit.py batch -p32 $case_bams \
      --normal $control_bams \
      --targets $intervals_list\
      --fasta $ref_fasta  \
      --output-reference my_reference.cnn --output-dir results \
      --scatter
  """
}
  

process TestProc {
    echo true

    input:
    file("*")
    file("*")
    
    """
    ls -al .
    """
}