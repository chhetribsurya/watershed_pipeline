#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PullSingularityImages {
    /*
     This process pulls specified Singularity images from Docker containers and
    stores them in a cache directory defined in the Nextflow config.
    */

    executor = 'local'
    debug true
    conda '../env_ymls/watershed_env.yml'
    publishDir "${params.singularity_cache_dir}", mode: 'copy'

    output:
        path("*.sif")        

    script:
    """
    #!/bin/bash
    set -e

    echo -e "\nChecking Singularity version and configuration..."
    singularity --version

    echo "Cache directory is set to: ${params.singularity_cache_dir}"
    echo -e "\nPulling Singularity images to cache..."
    #singularity pull --dir ${params.singularity_cache_dir} docker://ensemblorg/ensembl-vep
    singularity pull docker://bryancquach/peer:1.3
    singularity pull docker://ensemblorg/ensembl-vep
    echo "Images pulled successfully."

    echo -e "\nListing contents of cache directory post-pull..."
    ls
    """
}

/*
process build_vep_cache_test {
    debug true
    tag { "JOB VEPCACHE" }
    label 'vep_docker_env'
    conda '../../env_ymls/watershed_env.yml'
    publishDir "$params.cachedir_vep", mode: 'move'

    output:
        path("*"), emit: true, optional: true

    script:
        def outputdir = "$params.cachedir_vep/vep_cache" // Use an absolute path here
        //def outputdir = "cachedir/vep_cache"

        """
        #mkdir "$outputdir"
        echo -e "\nBuilding VEP cache GRCh38 ....\n"
        INSTALL.pl --CACHEDIR "$outputdir" -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE -PLUGINSDIR "$outputdir"/Plugins/

        #INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE -PLUGINSDIR ${params.cachedir_vep}/Plugins/
        """
}
*/


/*
process build_vep_cache_test2 {
    debug true
    tag { "JOB VEPCACHE" }
    label 'vep_docker_env'
    //publishDir "$params.cachedir_vep", mode: 'move'
    
    output:
        path("*"), emit: true, optional: true

    script:
        """
        mkdir -p /vepcache/vep_cache
        echo -e "\nBuilding VEP cache GRCh38 ....\n"
        INSTALL.pl --CACHEDIR /vepcache/vep_cache -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE -PLUGINSDIR /vepcache/vep_cache/Plugins/
        """
}
*/

process build_vep_cache {
    debug true
    tag { "JOB VEPCACHE" }
    label 'vep_docker_env'
    publishDir "$params.cachedir_vep", mode: 'move'
    
    output:
        path("*"), optional: true, emit: cache
        //path("*"), emit: true, optional: true

    script:
        """
        mkdir -p /vepcache
        echo -e "\nBuilding VEP cache GRCh38 ....\n"
        INSTALL.pl --CACHEDIR /vepcache -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE -PLUGINSDIR /vepcache/Plugins/
        """
}


process chromsplit_rv_file {
    debug true
    executor = 'local'
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/testvep", mode : 'copy'

    input:
        tuple val(pop), val(gene_rarevariant_file)

    output:
        path("splitvcf/${pop}/chromwise/*.txt"), emit: annotations
        path("splitvcf/${pop}/chromwise/*.vcf"), emit: vcfannotations

    script:
        def outputdir = "splitvcf/${pop}/chromwise"
        def script_dir = "${launchDir}/bin"
    
        """
        mkdir -p "$outputdir"
        echo "Splitting chromwise to VCF format: rarevariant file: $gene_rarevariant_file, output_dir: $outputdir, population : $pop"
        bash ${launchDir}/bin/step6_chromwise_rare_variants_to_vcf.sh "${pop}" "$gene_rarevariant_file" "$outputdir" "$script_dir"
        """
}


process call_vep_annotation_final {
    debug true
    tag { "JOB ${vcf_file}" }
    label 'vep_docker_env'
    publishDir path: "$params.out_dir/outputs/testvep", mode : 'copy'

    //publishDir path: "$params.out_dir/outputs/testvep", pattern: "gene*", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file)
        //path(vcf_file) 
        //path(cachedir_vep)

    output:
        path("vep/$pop/gene*"), optional: true, emit: all
        path("vep/$pop/gene*.vcf"), optional: true, emit: annotations
        path("vep/$pop/gene*.html"), optional: true, emit: html

    script:
        def outputdir = "vep/$pop"
        def basename = vcf_file.baseName.replace(".CADD", "")
        def outfile = "$outputdir/${basename}.vep.loftee.vcf"

        """
        mkdir -p "$outputdir"
        echo -e "\nAnnotation of VCFs with VEP ....\n"
        vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache /vepcache --plugin LoF,loftee_path:/vepcache/Plugins/loftee --dir_plugins /vepcache/Plugins/loftee
        """
}


//vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache ${params.cachedir_vep} --plugin LoF,loftee_path:${params.cachedir_vep}/Plugins/loftee --dir_plugins ${params.cachedir_vep}/Plugins/loftee



//-------------------------------------------------------------------------------------------------



workflow {

    PullSingularityImages()
    build_vep_cache()

    //Below handles the VEP chromwise run
    Channel
        .value("$params.cohort_name")
        .set { pops_annot_ch }

    Channel
        .fromPath('../data/gene*GLOBAL*txt')
        .set { gene_rarevar_ch } //annotation file

    combined_ch = pops_annot_ch
        .merge(gene_rarevar_ch)
        .map { pops_annot, gene_rarevar -> tuple(pops_annot, gene_rarevar)
    }

    combined_ch
        .view({ pops_annot, gene_rarevar  -> "Final Global annotation call for: ($pops_annot, $gene_rarevar)" }) 

    //Split files by chromosomes
    chromsplit_rv_file(combined_ch)

    vcffiles_ch = chromsplit_rv_file.out.vcfannotations

    /*vcffiles_ch
        .view({"Here's the VCF file list: $it"})*/

    vcffiles_ch.flatten()
        .take(params.num_chromosomes)
        .map { vcf -> tuple("$params.cohort_name", vcf) }
        .set { vcf_channel }

    //Run VEP Annotation chromosome wise as a separate process
    call_vep_annotation_final(vcf_channel)


}



//-------------------------------------------------------------------------------------------------


// IGNORE BELOW
/*process PullSingularityImages {
    executor = 'local'
    debug true
    //conda '../../env_ymls/watershed_env.yml'

    script:
    """
    #!/bin/bash
    set -e
    echo "Pulling Singularity images..."
    # Create the directory if it does not exist
    #mkdir -p ${params.singularity_cache_dir}
    #echo "Ensured cache directory exists."
    #ls -lah ${params.singularity_cache_dir}


    #singularity pull --dir ${params.singularity_cache_dir} docker://ensemblorg/ensembl-vep
    #singularity pull --dir ${params.singularity_cache_dir} docker://bryancquach/peer:1.3
    singularity pull docker://ensemblorg/ensembl-vep
    singularity pull docker://bryancquach/peer:1.3

    echo "Images pulled successfully."
    """
}
*/


/*process PullSingularityImages {
    debug true
    //
    // This process pulls specified Singularity images from Docker containers and
    // stores them in a cache directory defined in the Nextflow config.
    //

    conda '../../env_ymls/watershed_env.yml'

    script:
    """
    #!/bin/bash
    set -e
    echo "Pulling Singularity images..."

    # Pull each image using Singularity. This respects the cacheDir set in the config.
    #singularity --help
    singularity pull docker://ensemblorg/ensembl-vep
    singularity pull docker://bryancquach/peer:1.3

    echo "Images pulled successfully."
    """
}
*/

/*workflow {
    // Call the process to pull Singularity images
    PullSingularityImages()
}*/
