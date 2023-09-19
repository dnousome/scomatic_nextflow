#!/usr/bin/env nextflow
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )
STARfasta="/fdb/STAR_indices/2.7.10b/UCSC/hg38/ref.fa"
RNAediting="/data/nousomedr/SComatic/RNAediting/AllEditingSites.hg38.txt"
PON="/data/nousomedr/SComatic/PoNs/PoN.scRNAseq.hg38.tsv"
FILTEREDBED="/data/nousomedr/SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed"
scriptsplitcell="/data/nousomedr/SComatic/scripts/SplitBam/SplitBamCellTypes.py"
scriptbasecounter="/data/nousomedr/SComatic/scripts/BaseCellCounter/BaseCellCounter.py"
scriptmergecounts="/data/nousomedr/SComatic/scripts/MergeCounts/MergeBaseCellCounts.py"
scriptbasecall1="/data/nousomedr/SComatic/scripts/BaseCellCalling/BaseCellCalling.step1.py"
scriptbasecall2="/data/nousomedr/SComatic/scripts/BaseCellCalling/BaseCellCalling.step2.py"
outdir=file(params.output)


workflow {
    file_input=Channel.fromPath(params.file_input)
                        .splitCsv(header: true, sep: "\t",strip:true)   
                        .map{ row ->
                        tuple(row.sample,row.bam,file("${row.bam}.bai"),row.barcodes)
                        }       
    file_input.view()
  

   // calculate_barcodes(file_input)
    //barcode_read=Chanell.fromPath(file_input)
    //Count the number of bams that are separated by Barcodes
    scomatic_step1(file_input)
    
    scin1=scomatic_step1.out.transpose()
    .map{sample,bam,bai-> tuple(sample,
    bam.Name.split('\\.')[1],bam,bai)}
    scomatic_step2(scin1)
    scin2=scomatic_step2.out.groupTuple()
    scomatic_step3(scin2) | scomatic_step4.map{sample,step1,step2,pass -> tuple(sample,pass)} | annotate_pass

}



process calculate_barcodes{
    input:
    tuple val(sample), path(bam), path(bai), path(barcodes)
    
    output:
    tuple val(sample), path("${sample}.celltypes")
    
    shell:
    """
    awk -F '\\t' '{print \$2}' ${barcodes} | sort | uniq  >${sample}.celltypes
    
    #celltypes="\$(cat celltypes)"
    #echo "\${celltypes}"

    """
}


process scomatic_step1 {
    conda '/data/nousomedr/programs/mambaforge/envs/SComatic'
    
    publishDir("${outdir}/scomatic", mode: 'copy')

    input:
    tuple val(sample), path(bam), path(bai), path(barcodes)
    
    output:
        tuple val(sample), path("*.bam"), path("*.bai")
    shell :
    """
    python $scriptsplitcell --bam ${bam} --meta ${barcodes} --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --max_NH 1 

    """

    stub:
    """
    touch "${sample}.testcelltype1.bam"
    touch "${sample}.testcelltype1.bam.bai"
    touch "${sample}.testcelltype2.bam"
    touch "${sample}.testcelltype2.bam.bai"
    """
}

process scomatic_step2 {
    conda '/data/nousomedr/programs/mambaforge/envs/SComatic'
    
    publishDir("${outdir}/scomatic", mode: 'copy')

    input:
        tuple val(sample), val(celltype), path("${sample}.${celltype}.bam"), path("${sample}.${celltype}.bam.bai")
    
    output:
        tuple val(sample), val(celltype), 
        path("Step2/${sample}.${celltype}.tsv"), optional:true

    shell :
    """
    mkdir Step2

    ##STEP2
    python $scriptbasecounter --bam ${sample}.${celltype}.bam \
        --ref $STARfasta \
        --chrom all \
        --out_folder Step2 \
        --min_bq 30 \
        --nprocs 24

    """

    stub:
    """
    mkdir Step2

    touch "Step2/${sample}.${celltype}.tsv"
    """
}
    
process scomatic_step3 {
    conda '/data/nousomedr/programs/mambaforge/envs/SComatic'
    
    publishDir("${outdir}/scomatic", mode: 'copy')

    input:
    tuple val(sample), val(celltypes), path(tsvs)
    
    output:
        tuple val(sample), path("Step3/${sample}.BaseCellCounts.AllCellTypes.tsv")
  

    shell :
    """
    ##STEP3

    mkdir Step3

    python $scriptmergecounts --tsv_folder . \
    --outfile Step3/${sample}.BaseCellCounts.AllCellTypes.tsv
    
    """

    stub:
    """
    mkdir Step3

    touch "Step3/${sample}.BaseCellCounts.AllCellTypes.tsv"
    """
}

process scomatic_step4 {
    conda '/data/nousomedr/programs/mambaforge/envs/SComatic'
    
    publishDir("${outdir}/scomatic", mode: 'copy')

    input:
    tuple val(sample), path(allcelltsv)
    
    output:
    tuple val(sample),
    path("Step4/${sample}.calling.step1.tsv"), 
    path("Step4/${sample}.calling.step2.tsv"), 
    path("Step4/${sample}.calling.step2.pass.tsv")

    shell:
    """
    mkdir Step4

    ##STEP4
    python $scriptbasecall1 \
    --infile ${sample}.BaseCellCounts.AllCellTypes.tsv \
    --outfile Step4/${sample} \
    --ref $STARfasta

    #STEP4.2
    python $scriptbasecall2 \
    --infile Step4/${sample}.calling.step1.tsv \
    --outfile Step4/${sample} \
    --editing $RNAediting \
    --pon $PON
    
    ##Filtering through UCSC
    bedtools intersect -header -a Step4/${sample}.calling.step2.tsv -b $FILTEREDBED | awk '\$1 ~ /^#/ || \$6 == "PASS"' > Step4/${sample}.calling.step2.pass.tsv
    """

    stub:
    """
    mkdir Step4

    touch "Step4/${sample}.calling.step1.tsv" 
    touch "Step4/${sample}.calling.step2.tsv"
    touch "Step4/${sample}.calling.step2.pass.tsv"
    """
}


process annotate_pass {
    module = ['annovar/2020-06-08']

    publishDir("${outdir}/scomatic/annotated", mode: 'copy')

    input:
    tuple val(sample), path("Step4/${sample}.calling.step2.pass.tsv")

    output:
    tuple val(sample), path("${sample}_annotated")

    shell:

    '''
    grep -v '#' Step4/!{sample}.calling.step2.tsv  | awk -F '-' -v OFS='\t' '{print $1,$2,$3,$4,$5,$0}' > sample.variants.avinput 
    table_annovar.pl sample.variants.avinput !{ANNOVAR_DATA}/hg38 \
        --thread 4 \
        --buildver hg38 \
        --outfile !{sample}_annotated \
        --remove \
        --protocol refGene,clinvar_20230416,cosmic92_coding,avsnp150,dbnsfp42a \
        --operation g,f,f,f,f \
        --nastring '.' -polish -otherinfo

    '''
}

