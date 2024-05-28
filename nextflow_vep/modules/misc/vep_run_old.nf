#!/usr/bin/env nextflow

process VEP_RUN {

    tag { "JOB ${vcf_file}" }
    publishDir path: "${params.resultoutdir}/vep", pattern: "variant_effect*", mode : 'copy'

    //publishDir path: params.outdir, pattern: "variant_effect*", mode : 'move'
    //container "./singularity-images/vep.sif"   
 
    input:
        path(vcf_file)

    output:
       path("variant_effect*"), optional: true, emit: all
       path("variant_effect*.txt"), optional: true, emit: vcf
       path("variant_effect*.html"), optional: true, emit: html

    script:
        def outfile = "variant_effect_${vcf_file.baseName}.txt"


        if(params.skip_cache)
        """
            echo -e "\nSkipping VEP cache ....\n"
            echo -e "\nAnnotation of VCFs with VEP ....\n"
            vep -i ${vcf_file} --cache --dir_cache ${params.cachedir_vep} --format vcf --tab --fork 4 --cache --force --regulatory --variant_class --max_af --af_gnomadg --output_file ${outfile}

        """


        else if(params.mode == "help")
        """
            echo -e "\nPrinting VEP help message ....\n"
            vep --dir \$PWD/vep_data --help

        """


        else 
        """
            echo -e "\nBuilding VEP cache GRCh38 ....\n"
            INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE  -PLUGINSDIR ${params.cachedir_vep}/Plugins/
            echo -e "\nBuilding VEP cache GRCh37 ....\n"
            INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh37 --CONVERT --PLUGINS dbNSFP,CADD,Carol,NearestGene,CAPICE  -PLUGINSDIR ${params.cachedir_vep}/Plugins/
            echo -e "\nAnnotation of VCFs with VEP ....\n"
            vep -i ${vcf_file} --cache --dir_cache ${params.cachedir_vep} --format vcf --tab --fork 4 --cache --force --regulatory --variant_class --max_af --af_gnomadg --output_file ${outfile}

        """
}
