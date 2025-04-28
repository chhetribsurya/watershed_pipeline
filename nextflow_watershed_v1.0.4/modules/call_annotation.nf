#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process call_phylopscore_annotation {
    tag { "JOB ${pop}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(gene_rarevariant_file)

    output:
        path("ucsc/${pop}/*.bed"), emit: annotations
        val(pop), emit: pops_annotations

    script:
        def outputdir = "ucsc/${pop}"

        """
        python ${launchDir}/bin/step9.2_phylop100way_query.py --phylop "${params.phylopannot_file}" --rarevariant $gene_rarevariant_file --pop $pop --outdir $outputdir 
        """
}

process call_gencode_annotation_final {
    tag { "JOB ${pop}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(gene_rarevariant_file)

    output:
        path("gencode/${pop}/*gencode.txt"), emit: annotations
        val(pop), emit: pops_annotations

    script:
        def outputdir = "gencode/${pop}"
        def output_file = "$outputdir/gene-${pop}-rv.gencode.txt"
    
        """
        mkdir -p "$outputdir"
        python ${launchDir}/bin/step10_gencode_tss_tes_distance.py "${params.gencode_file}" "$gene_rarevariant_file" "$output_file"
        """
}


process chromsplit_rv_file {
    tag { "CADD-JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(gene_rarevariant_file)

    output:
        path("cadd/${pop}/chromwise/*.txt"), emit: annotations
        path("cadd/${pop}/chromwise/*.vcf"), emit: vcfannotations

    script:
        def outputdir = "cadd/${pop}/chromwise"
        def script_dir = "${launchDir}/bin"
    
        """
        mkdir -p "$outputdir"
        bash ${launchDir}/bin/step6_chromwise_rare_variants_to_vcf.sh "${pop}" "$gene_rarevariant_file" "$outputdir" "$script_dir"
        """
}

process call_cadd_annotation_final {
    debug true
    tag { "JOB $params.cohort_name" }
    label 'watershed_cadd_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        val(vcf)

    output:
        path("cadd/${params.cohort_name}/*tsv.gz"), emit: annotations
        val("$params.cohort_name"), emit: pops_annotations


    script:
        def pop = "$params.cohort_name"
        def outputdir = "cadd/${pop}"
        def vcf_name = "$vcf.baseName"
        def chrom = vcf_name.split('\\.').find { it.startsWith('chr') }
        def output_file = "$outputdir/gene-${pop}-rv.CADD.${chrom}.tsv.gz"


        def cadd_script_files = "$params.cadd_script_files"
        def corenum = 8
        def genome = "GRCh38"
        def inputvcf_file = "$vcf"
    
        """
        mkdir -p "$outputdir"
        bash ${cadd_script_files}/CADD.sh -a -g ${genome} -c ${corenum} -o ${output_file} ${inputvcf_file}
        """
}

process call_cadd_annotation {
    tag { "JOB ${pop}" }
    label 'watershed_cadd_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(vcf)

    output:
        path("cadd/${pop}/*tsv.gz"), emit: annotations
        val(pop), emit: pops_annotations
        val("$params.out_dir/outputs/2_annotation/cadd/${pop}"),  emit: cadd_dir 


    script:
        def outputdir = "cadd/${pop}"
        def vcf_name = "$vcf.baseName"
        def chrom = vcf_name.split('\\.').find { it.startsWith('chr') }
        def output_file = "$outputdir/gene-${pop}-rv.CADD.${chrom}.tsv.gz"

        def cadd_script_files = "$params.cadd_script_files"
        def corenum = 8
        def genome = "GRCh38"
        def inputvcf_file = "$vcf"
    
        """
        mkdir -p "$outputdir"
        bash ${cadd_script_files}/CADD.sh -a -g ${genome} -c ${corenum} -o ${output_file} ${inputvcf_file}
        """
}

/*process call_vep_annotation_test {

    tag { "JOB ${vcf_file}" }
    publishDir path: "$params.out_dir/outputs/2_annotation", pattern: "gene*", mode : 'copy'
 
    input:
        path(vcf_file)

    output:
       path("vep/${pop}/gene*"), optional: true, emit: all
       path("vep/${pop}/gene*.vcf"), optional: true, emit: annotations
       path("vep/${pop}/gene*.html"), optional: true, emit: html

    script:
        def outputdir = "vep/${pop}"
        def newbasename = vcf_file.baseName.replace(".CADD", "")
        def outfile = "$outputdir/${newbasename}.vep.loftee.vcf"

        if(params.skip_cache)
        """
            mkdir -p "$outputdir"

            echo -e "\nSkipping VEP cache ....\n"
            echo -e "\nAnnotation of VCFs with VEP ....\n"
            vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache ${params.cachedir_vep} --plugin LoF,loftee_path:${params.cachedir_vep}/Plugins/loftee --dir_plugins ${params.cachedir_vep}/Plugins/loftee

        """

        else if(params.mode == "help")
        """
            echo -e "\nPrinting VEP help message ....\n"
            vep --help
            //vep --dir \$PWD/vep_data --help

        """

        else 
        """
            mkdir -p "$outputdir"

            echo -e "\nBuilding VEP cache GRCh38 ....\n"
            INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE  -PLUGINSDIR ${params.cachedir_vep}/Plugins/

            //echo -e "\nBuilding VEP cache GRCh37 ....\n"
            //INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh37 --CONVERT --PLUGINS dbNSFP,CADD,Carol,NearestGene,CAPICE  -PLUGINSDIR ${params.cachedir_vep}/Plugins/

            echo -e "\nAnnotation of VCFs with VEP ....\n"
            vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache ${params.cachedir_vep} --plugin LoF,loftee_path:${params.cachedir_vep}/Plugins/loftee --dir_plugins ${params.cachedir_vep}/Plugins/loftee

            //vep -i ${vcf_file} --cache --dir_cache ${params.cachedir_vep} --format vcf --tab --fork 4 --cache --force --regulatory --variant_class --max_af --af_gnomadg --output_file ${outfile}
        """
}
*/

process build_vep_cache {
    tag { "JOB VEPCACHE" }
    label 'watershed_env'
    label 'vep_docker_env'
    publishDir "$params.cachedir_vep", mode: 'copy'

    output:
        path("cachedir/vep_cache/*"), emit: cache_files

    script:
        def outputdir = "cachedir/vep_cache"

        """
        mkdir "$outputdir"
        INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE -PLUGINSDIR ${params.cachedir_vep}/Plugins/
        """
}


process call_vep_annotation_final {
    tag { "JOB $params.cohort_name" }
    label 'vep_docker_env'
    publishDir path: "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file)

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
        vep -i ${vcf_file} --format vcf --output_file ${outfile} --vcf --cache --regulatory --fork 4 --cache --dir_cache /vepcache --plugin LoF,loftee_path:/vepcache/Plugins/loftee --dir_plugins /vepcache/Plugins/loftee
        """
}


process prepare_vep_parse {
    tag { "JOB $params.cohort_name" }
    label 'watershed_vepenv'
    publishDir path: "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        each vcf_file
        val(pop)

    output:
        path("vep/$pop/gene-${pop}-rv*vep.loftee*vcf.gz"), emit: vep_vcf_annotations
        path("vep/$pop/gene-${pop}-rv*vep.loftee*vcf.gz*"), emit: all_annotations

    script:
        def outputdir = "vep/$pop"
        def outfile = "$outputdir/gene-${pop}-rv.vep.loftee.snv.tsv"

        def vcf = vcf_file.baseName
        def gz_file = "$outputdir/${vcf}.vcf.gz"

        """
        mkdir -p $outputdir
        bgzip -c $vcf_file > $gz_file
        tabix -p vcf $gz_file
        ls -al $outputdir
        """
}

process parse_vep_annotation {
    tag { "JOB $params.cohort_name" }
    label 'watershed_vepenv'
    publishDir path: "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file_list)

    output:
        path("vep/$pop/gene-${pop}-rv.vep.loftee.snv.tsv"), emit: parsed_annotations

    script:
        def outputdir = "vep/$pop"
        def outfile = "$outputdir/gene-${pop}-rv.vep.loftee.snv.tsv"
        def vcf_files_string = vcf_file_list.join(',')

        """
        mkdir -p $outputdir
        python ${launchDir}/bin/step11.1_parse_vep_loftee_chromwise-final.py --anno_vcf_list $vcf_files_string --outputfile $outfile
        """
}

