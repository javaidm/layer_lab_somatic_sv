import static nextflow.Nextflow.file
// import static nextflow.Nextflow.file

import nextflow.Channel

class LLabUtils {

static def extractSamples(tsvFile) {
  // Channeling the TSV file containing FASTQ or BAM
  // Format is: "subject gender status sample lane fastq1 fastq2"
  // or: "subject gender status sample lane bam"

  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    def parent_dir = file(tsvFile).parent
    def sample  = row[0].trim()
    def trimmed_r1 = "${row[1]}".trim()
    def trimmed_r2 = "${row[2]}".trim()
    def path1 = "${parent_dir}/${trimmed_r1}"
    def path2 = "${parent_dir}/${trimmed_r2}"
    // println("path1:$path1")
    // println("path2:$path2")
    def file1      = this.returnFile(path1)
    def file2      = this.returnFile(path2)
    [sample,[ file1, file2]]
  }
  
} // end of method extractSample()

static def extractControls(tsvFile) {
  def normals = []
  def allLines = tsvFile.readLines()
  for (line in allLines){
    def trimmed = line.trim()
    def columns = trimmed.split()
    
    // pick the first column value for normals
    if (columns.size() > 0 ){
      normals.add(columns[0])
    }    
  }
  return normals
} // end of method extractControls()

static def sampleInList(sample_file_name, list_of_samples){
  def retVal = false
  list_of_samples.each{
      // println("checking $it")
      if (sample_file_name ==~ /.*$it.*/ )
        retVal = true
     }
     return retVal
}


static def extractCases(tsvFile) {
  // assuming the first column in the control, and the rest are the cases
  def cases = []
  def allLines = tsvFile.readLines()
  for (line in allLines){
    def trimmed = line.trim()
    def columns = trimmed.split()
    columns.eachWithIndex { it, i ->
      if (i != 0){
        // println("i: $i")
        cases.add(it)
      }
    }
  } // end for
  return cases
} // end of method extractCases()

static def extractSampleAndConditionFiles(crams_dir){
  def process_list=[]
  def p_crams = file("$crams_dir/*P.cram")
  
  for( def f : p_crams ) {
    //extract sample name
    def index = f.baseName.indexOf('P.cram')
    def normal_sample_name = f.baseName
    def sample = f.baseName[0..(index-1)]
    //get sample name for the C, G, and GC file
    def file_list = file("$crams_dir/${sample}C_TD*.cram")    
    process_list.add(
      [sample, normal_sample_name, 'C_P', [file_list[0], f]]
      )
    
    file_list = file("$crams_dir/${sample}G_TD*.cram")    
    process_list.add(
      [sample, normal_sample_name, 'G_P',[file_list[0], f]]
      )

    file_list = file("$crams_dir/${sample}GC_TD*.cram")    
    process_list.add(
      [sample, normal_sample_name, 'GC_P',[file_list[0], f]]
      )
  }
  // return Channel.from(process_list)
  return process_list
}
// static def extractSampleAndConditionFiles(crams_dir, file_type){
//   def process_list=[]
//   def p_crams = file("$crams_dir/*P_TD*.cram")
  
//   for( def f : p_crams ) {
//     //extract sample name
//     def index = f.baseName.indexOf('P_TD')
//     def sample = f.baseName[0..(index-1)]
//     //get sample name for the C, G, and GC file
//     def file_list = file("$crams_dir/${sample}${file_type}_TD*.cram")    
//     process_list.add([sample, [file_list[0], f]])
//   }
//   // return Channel.from(process_list)
//   return process_list
// }

static def runSanityChecks(dirPath){
  // this.checkIfBothReadsExist(dirPath)
}




// Layer Lab ascii art
  static def layerLabAscii() {
    // def ascii_str = 
    println \
    '''
    | |                         | |         | |    
    | |     __ _ _   _  ___ _ __| |     __ _| |__ 
    | |    / _` | | | |/ _ \\ '__| |    / _` | '_ \\ 
    | |___| (_| | |_| |  __/ |  | |___| (_| | |_) |
    |______\\__,_|\\__, |\\___|_|  |______\\__,_|_.__/ 
                __/ |                            
               |___/  
    '''
    // println ascii_str
  }

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  static def getChrmList(){
    def chrs = (1..22).collect()
    chrs.addAll(['X', 'Y', 'MT'])
    return chrs
  }

  static def getChrmListHg38(){
    def chrs = (1..22).collect()
    chrs.addAll(['X', 'Y'])
    return chrs.collect {"chr$it"}
  }
  
  // // Return element in list of allowed params
  static def checkParams(it) {
    return it in [
      'contact-mail',
      'contactMail',
      'container',
      'dbsnp',
      'gff3',
      'help',
      'known-indels',
      'knownIndels',
      'results-dir',
      'resultsDir',
      'publish-dir-mode',
      'publishDirMode',
      'ref-fasta',
      'refFasta',
      'run-name',
      'runName',
      'sample-dir',
      'sampleDir',
      'sample-tsv',
      'sampleTsv',
      'tag',
      'verbose',
      'version'
      ]
  }

  // // Check file extension
  static def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
  }

  static def exists (path){
    return file(path).exists()
  }
  // // Return file if it exists
  static def returnFile(it) {
    // if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    if (!this.exists(it)) exit 1, "Missing file: ${it}, see --help for more information"
    return file(it)
  }

  // // Return status [0,1]
  // // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
  }
}
