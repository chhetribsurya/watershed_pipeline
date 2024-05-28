#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Input parameters: read pairs
 * Params are stored in the params.config file
 */

version                 = "1.0"

// this prevents a warning of undefined parameter
//params.help            = false
params.genome            = false
params.reads             = false


println """\

         VEP ANNOTATION - N F   P I P E L I N E  ~  Version ${version}
         ==========================================================================================
         #genome      : ${params.genome}
         #reads       : ${params.reads}
         inputvcfs    : ${params.inputvcf}
         outputdir    : ${params.outdir}
         cachedir     : ${params.cachedir}
         resultdir    : ${params.resultoutdir}

         """
//         .stripIndent()



/*
===================================================================================================
    Include Modules
===================================================================================================
*/

include { VEP_RUN_CHROMWISE }                                from "./modules/vep_run_chromwise" addParams(OUTPUT: "${params.outdir}/vep_out")

/*
include { VEP_RUN    }                                       from "./modules/vep_run" addParams(OUTPUT: "${params.outdir}/vep_out")
*/


/*
==================================================================================================
    Create Channels
==================================================================================================
*/

Channel
    .fromPath(params.inputvcf, checkIfExists: true)
    .set({vcf_data})



/*
==================================================================================================
    WORKFLOW - VEP Annotation for Watershed
==================================================================================================
*/

// Define parameters
params.vcfDir = "./chromwise"

// Define the input channel
input_ch = Channel.fromPath("${params.vcfDir}/gene-*-rv.chr*.vcf")

input_ch.view({ "Found file: $it" })

workflow{

    //VEP_RUN(vcf_data)
    VEP_RUN_CHROMWISE(input_ch)

    //results_vep_ch.view({"Here's the VCF process output file: $it"})
    //VEP_RUN.out.view({"\n" + "Here's the VCF process output file -- Case1: $it" + "\n"}) //outputs everything

    //specific emit based output
    VEP_RUN_CHROMWISE.out.all.view({"\n" + "Here's the all VCF process output file -- Case2: $it" + "\n"}) //outputs all (everything)
    VEP_RUN_CHROMWISE.out.vcf.view({"\n" + "Here's the only VCF process output file -- Case3: $it" + "\n"}) //outputs vcf files
    VEP_RUN_CHROMWISE.out.html.view({"\n" + "Here's the only HTML process output file -- Case4: $it" + "\n"}) //outputs html files

}


/*
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
*/

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}


/*
========================================================================================
    THE END
========================================================================================
*/
