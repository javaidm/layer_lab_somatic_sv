# Just a sample makefile to show one example run of the pipeline
all:
	nextflow run main.nf \
		--sample-tsv=/Shares/layer_shared/projects/cancer_center_tiny/samples.tsv \
		--conditions-tsv /Shares/layer_shared/projects/cancer_center_tiny/conditions.tsv \
		--run-name='cancer_center_tiny' -profile fiji -resume

yaml:
	nextflow run main.nf -params-file params.yaml -profile fiji -resume