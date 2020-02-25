#!/usr/bin/env nextflow
nextflow.preview.dsl=2
initParamsToDefaults()
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
  gnomad_genome = params.gnomadGenome?: "$fiji_gnomad_genome"
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
    inputFiles = LLabUtils.extractSample(tsvFile)
}else{
  // Define channel for reading file pairs
  Channel
      .fromFilePairs( params.reads, size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}"}
      .set {inputFiles}
}
ch_conditions = LLabUtils.extractConditions(file(params.conditionsTsv))
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

Channel.from(LLabUtils.getChrmList())
.set{chrmList}
// CRAMS_DIR_ABS_PATH = "${outDir}/align"

workflow{
    // GATK portion
    PreprocessIntervals()
    CollectReadCounts(ch_crams,
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
    
  //  SavvyCNV related Summaries
   GenSavvyCNVCoverageSummary(ch_bams)
   RunSavvyCNV(GenSavvyCNVCoverageSummary.out.collect())
  
  // CNVKIT related portion
  RunCNVKitBatchMode(BAMS_DIR)
  //  PrepareBaitedTargets_cnvkit()
  //  GenAccessibleRegions_cnvkit()
  //  DeriveAntitargetRegions_cnvkit(PrepareBaitedTargets_cnvkit.out,
  //                           GenAccessibleRegions_cnvkit.out)
  //  RunAutoBin_cnvkit(ch_all_bams.collect(),
  //                           PrepareBaitedTargets_cnvkit.out,
  //                           GenAccessibleRegions_cnvkit.out)
  //  GetCoverage_cnvkit(ch_bams,
  //                           PrepareBaitedTargets_cnvkit.out,
  //                           DeriveAntitargetRegions_cnvkit.out)
  // GetPooledRef_cnvkit(GetCoverage_cnvkit.out[0].collect(),
  //                     GetCoverage_cnvkit.out[1].collect())
  // GetCnr_cnvkit(GetCoverage_cnvkit.out[0],
  //                     GetCoverage_cnvkit.out[1],
  //                     GetPooledRef_cnvkit.out)
  // GetSegments_cnvkit(GetCnr_cnvkit.out)
  // CallSegments_cnvkit(GetSegments_cnvkit.out)
  // PlotScatter_cnvkit(GetCnr_cnvkit.out,
  //                   GetSegments_cnvkit.out)

  // // PlotIdeogram_cnvkit(GetCnr_cnvkit.out,
  // //                   GetSegments_cnvkit.out)
                    
  // PlotHeatmap_cnvkit(GetSegments_cnvkit.out.collect())

} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

def initParamsToDefaults(){
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
  runStr = "nextflow run javaidm/layer_lab_dna_seq_vc "
  log.info "    Usage:"
  log.info "       $runStr --sample <file.tsv>"
  log.info "       $runStr --reads <shell glob pointing to the reads>"
  log.info "    --sample-tsv <file.tsv>"
  log.info "       Specify a TSV file containing sample_id, path_to_read1, path_to_read2."
  log.info "    --reads <shell glob>"
  log.info "       Specify a shell glob pointing to the reads."
  log.info "    --results-dir <directory>"
  log.info "       Specify a directory to hold the analysis results."
  log.info "    --run-name"
  log.info "       Specify a run name, results will be stored under results-dir/run-name"
  log.info "    --manifest [optional <manifest.csv>]"
  log.info "       Specify an optional manifest file."  
  log.info "    --dbsnp [optional if already specified as part of a Nextflow Profile]"
  log.info "    --ensembl-gene-annotation [optional if already specified as part of a Nextflow Profile]"
  log.info "    --known-indels [optional if already specified as part of a Nextflow Profile]"
  log.info "    --ref-fasta [optional if already specified as part of a Nextflow Profile]"
  log.info "       Human Genome Reference File."  
  log.info "    --onlyQC"
  log.info "       Run only QC tools and gather reports."
  log.info "    --help"
  log.info "       you're reading it."
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow."
}

/* Processes */
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
    // set val(sample), file (cram_index_pair)
    tuple val(sample),file(cram_index_pair)
    file(preprocessed_intervals)

    output:
    file(out_file)

    script:
    cram = cram_index_pair[0]
    cram_index = cram_index_pair[1]
    out_file = "${cram.simpleName}.counts.hdf5"
    
    
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
        if (it =~ /.*P\.counts\.hdf5/){
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
    cache false
    publishDir "$publish_dir", mode: 'copy', overwrite: true
    
    input:
    set val(sample), file(bam_index_pair)
    
    output:
    file(out_file)
    // file(publish_dir)

    script:
    // sample = bam.simpleName
    bam = bam_index_pair[0]
    publish_dir = "${OUT_DIR}/misc/SavvyCNV_coverage_summaries"
    out_file = "${sample}.coverageBinner"
    
    """
    java -Xmx1g CoverageBinner $bam >$out_file 
    """
}

process RunSavvyCNV {
    echo true
    cache false
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
  cache false
  publishDir "$publish_dir", mode: 'copy', overwrite: true
  input:
  path(bams_dir)
  output:
  file('results')
  script:
  publish_dir = "${OUT_DIR}/misc/CNVKit_calls"
  """
  cnvkit.py batch -p32 ${bams_dir}/*{C,G}.bam \
      --normal ${bams_dir}/*P.bam \
      --targets $intervals_list\
      --fasta $ref_fasta  \
      --output-reference my_reference.cnn --output-dir results \
      --scatter
  """
}
  
