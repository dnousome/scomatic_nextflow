# scomatic_nextflow
Nextflow implementation of the [SComatic](https://github.com/cortes-ciriano-lab/SComatic) somatic variant calling pipeline for scRNA data on the NIH HPC Biowulf system.

Files required for running are 1) The barcodes per sample for each celltype by sample and 2) a tab delimited sample sheet indicating the sample, bam/bai file location, and the barcode file location

barcodes.tsv

    Index   Cell_type
    AAATGGAAGGGACAGG        epithelial
    AACAAGAAGGTTGCCC        epithelial

samplesheet.tsv

    sample  bam     barcodes
    SAMN18434384    /data/darryl/SAMN18434384Aligned.sortedByCoord.out.bam      /data/darryl/SComatic/SAMN18434384_Barcodes.tsv

Code

    module load nextflow
    nextflow run scomatic.nf -c nextflow.config -profile biowulf -resume --file_input samplesheet.tsv --output sc_out1 
