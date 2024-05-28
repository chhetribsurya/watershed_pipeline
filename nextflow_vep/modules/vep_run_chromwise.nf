#!/usr/bin/env nextflow

process VEP_RUN_CHROMWISE {

    tag { "JOB ${vcf_file}" }
    publishDir path: "${params.resultoutdir}/vep", pattern: "gene*", mode : 'copy'

    //publishDir path: params.outdir, pattern: "variant_effect*", mode : 'move'
    //container "./singularity-images/vep.sif"   
 
    input:
        path(vcf_file)

    output:
       path("gene*"), optional: true, emit: all
       path("gene*.vcf"), optional: true, emit: vcf
       path("gene*.html"), optional: true, emit: html

    script:
        def newbasename = vcf_file.baseName.replace(".CADD", "")
        def outfile = "${newbasename}.vep.loftee.vcf"


        if(params.skip_cache)
        """
            echo -e "\nSkipping VEP cache ....\n"
            echo -e "\nAnnotation of VCFs with VEP ....\n"
            #vep -i ${vcf_file} --cache --dir_cache ${params.cachedir_vep} --format vcf --tab --fork 4 --cache --force --regulatory --variant_class --max_af --af_gnomadg --output_file ${outfile}
            #vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache ${params.cachedir_vep} --plugin LoF,loftee_path:/scratch16/abattle4/surya/datasets/watershed_rarevars/project_final/project-watershed/cachedir/vep_cache/Plugins/loftee --dir_plugins /scratch16/abattle4/surya/datasets/watershed_rarevars/project_final/project-watershed/cachedir/vep_cache/Plugins/loftee

            vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache ${params.cachedir_vep} --plugin LoF,loftee_path:${params.cachedir_vep}/Plugins/loftee --dir_plugins ${params.cachedir_vep}/Plugins/loftee

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
