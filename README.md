# Somatic Structural Variants Analysis Pipeline
 `PLEASE NOTE THAT THIS IS VERY MUCH A WORK IN PROGRESS AND THE CONFIGURATION AND POSSIBLY THE CODE MIGHT NEED TWEAKING BEFORE RUNNING THE PIPELINE`

The pipeline is based on [nextflow](https://www.nextflow.io) and uses the following SV callers:
1. [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs-) (Mutect2 based workflow with small panel of normals)
2. [SavvyCNV](https://github.com/rdemolgen/SavvySuite) (based on off-read signals)
3. [CNVKit](https://cnvkit.readthedocs.io/en/stable/)

## Pipeline Architecture
Below are the major compoents of the pipeline
- Nextflow is the main glue that wraps the individual commands and shell scripts (bwa mem for example), and orchestrate the workflow
- Configuration files (under the `conf` directory) carries the information that nextflow uses to work out which *executor* to use to run the pipeline
- Herlper utitlites (under `lib`) are some herlper funcions written to aid the main pipeline
- Individual scripts wrapped inside `nextflow processes` to carry out the actual tasks such as the *alighnemnt* or *marking duplicate reads*, or running a *structural variant caller*

## How to run the pipeline
For a run on your machine (or any other infrastructure such as an HPC, or AWS), you need to add corrresponding configuration under the `conf`. The *configuration* defines a *profile* in Nextflow lingo and needs to be passed at the commandline when running the pipeline. See the `Makefile` in the top level directory for an example run. At the command prompt (while running thorugh Nextflow), you will need to pass a *samples.tsv* (tab delimited) carrying atleast three columns. The first column specifying the *sample name*, and the next two specifying the paths to the *first* and the *second reads*.


