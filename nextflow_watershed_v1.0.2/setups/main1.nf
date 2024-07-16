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


//-------------------------------------------------------------------------------------------------



workflow {
    PullSingularityImages()
    build_vep_cache()
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

