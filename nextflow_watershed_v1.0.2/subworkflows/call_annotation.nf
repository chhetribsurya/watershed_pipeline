#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*include { process_data; call_phylopscore_annotation } from "../modules/call_annotation_process.nf" params(params)*/

process call_phylopscore_annotation {
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_python_env'
    publishDir "$params.out_dir/outputs/2_annotation", mode : 'copy'

    input:
        tuple val(pop), val(gene_rarevariant_file)

    output:
        path("ucsc/${pop}/*.bed"), emit: annotations
        val(pop), emit: pops_annotations

    script:
        /*def outputdir = "ucsc/${pop}"*/
        def outputdir = "ucsc/${pop}"

        """
        echo "phylop file: ${params.phylopannot_file}, rarevariant file: $gene_rarevariant_file, population : $pop"
        echo "querying phylop 100 way scores for population: $pop"
        python ${launchDir}/bin/step9.2_phylop100way_query.py --phylop "${params.phylopannot_file}" --rarevariant $gene_rarevariant_file --pop $pop --outdir $outputdir 
        """
}

process call_gencode_annotation_final {
    debug true
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
        echo "gencode file: ${params.gencode_file}, rarevariant file: $gene_rarevariant_file, output_file: $output_file, population : $pop"
        python ${launchDir}/bin/step10_gencode_tss_tes_distance.py "${params.gencode_file}" "$gene_rarevariant_file" "$output_file"
        """
}


process chromsplit_rv_file {
    debug true
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
        echo "Splitting chromwise to CADD VCF format: rarevariant file: $gene_rarevariant_file, output_dir: $outputdir, population : $pop"
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
        //tuple val(pop), val(vcf)

    output:
        path("cadd/${params.cohort_name}/*tsv.gz"), emit: annotations
        val("$params.cohort_name"), emit: pops_annotations
        //path("cadd/${pop}/*tsv.gz"), emit: annotations
        //val(pop), emit: pops_annotations


    script:
        def pop = "$params.cohort_name"
        def outputdir = "cadd/${pop}"
        def vcf_name = "$vcf.baseName"
        def chrom = vcf_name.split('\\.').find { it.startsWith('chr') }
        def output_file = "$outputdir/gene-${pop}-rv.CADD.${chrom}.tsv.gz"

        /*def output_file = "${outdir}/gene-${pop}-rv.CADD.${chrom}.tsv.gz"*/

        def cadd_script_files = "$params.cadd_script_files"
        def corenum = 8
        def genome = "GRCh38"
        def inputvcf_file = "$vcf"
    
        """
        mkdir -p "$outputdir"
        echo "CADD Command line: bash ${cadd_script_files}/CADD.sh -a -g ${genome} -c ${corenum} -o ${output_file} ${inputvcf_file}"
        bash ${cadd_script_files}/CADD.sh -a -g ${genome} -c ${corenum} -o ${output_file} ${inputvcf_file}
        """
}

process call_cadd_annotation {
    debug true
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
        echo "CADD Command line: bash ${cadd_script_files}/CADD.sh -a -g ${genome} -c ${corenum} -o ${output_file} ${inputvcf_file}"
        bash ${cadd_script_files}/CADD.sh -a -g ${genome} -c ${corenum} -o ${output_file} ${inputvcf_file}
        """
}

/*process call_vep_annotation_test {

    tag { "JOB ${vcf_file}" }
    publishDir path: "$params.out_dir/outputs/2_annotation", pattern: "gene*", mode : 'copy'
    //container "./singularity-images/vep.sif"   
 
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
    debug true
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
        echo -e "\nBuilding VEP cache GRCh38 ....\n"
        INSTALL.pl --CACHEDIR ${params.cachedir_vep} -a cfp -s homo_sapiens -y GRCh38 --CONVERT --PLUGINS dbNSFP,CADD,G2P,TSSDistance,FunMotifs,GWAS,LOEUF,Gwava,LoF,LoFtool,pLI,Conservation,LD,SpliceAI,SpliceRegion,REVEL,LOFTEE -PLUGINSDIR ${params.cachedir_vep}/Plugins/
        """
}


process call_vep_annotation_final {
    debug true
    //tag { "JOB ${vcf_file}" }
    tag { "JOB $params.cohort_name" }
    label 'vep_docker_env'
    publishDir path: "$params.out_dir/outputs/2_annotation", mode : 'copy'

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


process prepare_vep_parse {
    debug true
    //tag { "JOB ${vcf_file}" }
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
        echo "Received VCF file for ${pop}: $vcf_file"  // Debug statement to show all VCF files received
        bgzip -c $vcf_file > $gz_file
        tabix -p vcf $gz_file
        echo "List of processed files"
        ls -al $outputdir
        """
}

process parse_vep_annotation {
    debug true
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
        //IFS=',' read -r -a array <<< "$vcf_files_string" //unix
        //echo "${array[@]}"

        """
        mkdir -p $outputdir
        echo "Received Input file: $vcf_files_string"
        python ${launchDir}/bin/step11.1_parse_vep_loftee_chromwise-final.py --anno_vcf_list $vcf_files_string --outputfile $outfile
        """
}


workflow CALL_ANNOTATION {
    take:
    pops_annot_ch
    gene_rarevar_ch
    gene_rarevar_ch2


    main:
    if (params.analysis == "GLOBAL") {
        
        /*Channel.value("$params.cohort_name").set { pops_annot_ch }
        Channel.fromPath('./data/gene*GLOBAL*bed').set { gene_rarevar_ch } //annotation file
        Channel.fromPath('./data/gene*GLOBAL*txt').set { gene_rarevar_ch2 } //annotation file
        */

        combined_ch = pops_annot_ch
            .merge(gene_rarevar_ch)
            .map { pops_annot, gene_rarevar -> tuple(pops_annot, gene_rarevar)
        }

        combined_ch2 = pops_annot_ch
            .merge(gene_rarevar_ch2)
            .map { pops_annot, gene_rarevar -> tuple(pops_annot, gene_rarevar)
        }

        combined_ch.view({ pops_annot, gene_rarevar  -> "Final Global annotation call for: ($pops_annot, $gene_rarevar)" }) 
        
        // Phylop annotation
        call_phylopscore_annotation(combined_ch)

        // Gencode annotation
        call_gencode_annotation_final(combined_ch2)

        // Split VCF files by chrom
        chromsplit_rv_file(combined_ch2)

        // Split files output
        vcffiles_ch = chromsplit_rv_file.out.vcfannotations
        vcffiles_ch.view({"Here's the VCF file example: $it[0]"})

        vcffiles_ch.flatten()
            .take(params.num_chromosomes)
            .map { vcf -> tuple("$params.cohort_name", vcf) }
            .set { vcf_channel }

        /*Run CADD Annotation chromosome wise as a separate process*/
        call_cadd_annotation(vcf_channel)
        
        /*Run VEP Annotation chromosome wise as a separate process*/
        call_vep_annotation_final(vcf_channel)

        //Prepare VEP parse of VCFs
        Channel.value("$params.cohort_name").set { pop_channel}    
        vep_vcf_ch = call_vep_annotation_final.out.annotations
        vep_vcf_ch
            .collect()
            .set { vcf_list_channel }

        prepare_vep_parse(vcf_list_channel, pop_channel)

        //Parse VEP annotations
        vcf_newlist_channel = prepare_vep_parse.out.vep_vcf_annotations
        vcf_newlist_channel
            .collect()
            .map { files -> tuple("$params.cohort_name", files) }
            .set { combined_vcf_channel }

        parse_vep_annotation(combined_vcf_channel)


       /*------------------------------------------------------------------*/
       /* if (!params.skip_cache) {
            println "BUILDING VEP CACHE ..." 
            build_vep_cache()
            build_vep_cache
                .out
                .cache_files
                .set { cache_dir_ch }
        } else {
            println "SKIPPING VEP CACHE ..."
            Channel
                .fromPath(params.cachedir_vep)
                .set { cache_dir_ch }
        }*/

        /*Channel
            .fromPath('path/to/vcfs/*.vcf')  // Adjust this path to your VCF files
            .flatten()
            .set { vcf_channel }*/

        //vcffiles_ch.flatten().set { vcf_channel }

        /*vcf_channel
           .combine(cache_dir_ch)
           .view { combined -> "VEP CACHE DIR and VCF FILE COMBINED TUPLE: ${combined}" }
           .set { vep_input_ch }*/

        // Map to create tuples instead of combine
        /*vcf_channel
            .map { vcf -> tuple(vcf, cache_dir_ch) }
            .view { combined -> "VEP CACHE DIR and VCF FILE COMBINED TUPLE: ${combined}" }
            .set { vep_input_ch }*/

        //call_vep_annotation_final(vep_input_ch)
        /*call_vep_annotation_final(vcffiles_ch.flatten(), cache_dir_ch )*/

       /*------------------------------------------------------------------*/

    } else if (params.analysis == "STRATIFIED") {

        // Joining channels based on the key
        combined_ch = residual_expr_keyed
            .join(gene_rarevar_keyed)
            .map { key, residual_expr, gene_rarevar ->
                tuple(key, residual_expr, gene_rarevar)
            }
            .filter { key, residual_expr, gene_rarevar -> key == "EUR" || key == "AFR" }
            .view({ key, residual_expr, gene_rarevar -> "Final generare-variant combined tuple: ($key, $residual_expr, $gene_rarevar)" })
        
        // Call expression outliers using the synchronized combined channel
        call_expr_outliers(combined_ch)
        call_expr_outliers.out.outliers.view({"Here's the outlier file list: $it"})
        call_expr_outliers.out.pops_outliers.view({"Here's the outlier pop: $it"})

    }

    emit:
    phylop_annot_file  = call_phylopscore_annotation.out.annotations
    gencode_annot_file  = call_gencode_annotation_final.out.annotations
    chromsplit_file  = chromsplit_rv_file.out.annotations
    chromsplit_vcf_file  = chromsplit_rv_file.out.vcfannotations
    cadd_annot_file  = call_cadd_annotation.out.annotations
    vep_annot_file  = call_vep_annotation_final.out.annotations
    vep_annot_bgzipped = prepare_vep_parse.out.vep_vcf_annotations
    vep_annot_parsed_file = parse_vep_annotation.out.parsed_annotations
    cadd_annot_dir  = call_cadd_annotation.out.cadd_dir
    //vep_cache_file = build_vep_cache.out.cache_files
}


