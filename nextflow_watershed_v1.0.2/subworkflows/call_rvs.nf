#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*include { process_data; call_phylopscore_annotation } from "../modules/call_annotation_process.nf" params(params)*/

process splitvcfs_by_chrom {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file)

    output:
        path "splitvcfs/${pop}/*.vcf.gz", emit: splitvcfs
        path "splitvcfs/${pop}/*.vcf.gz.tbi", emit: indexes

    script:
        def outputdir = "splitvcfs/${pop}"
        //def outfile = "${outputdir}/gene-${pop}-rv.N2pairs.tsv"

        """
        #!/usr/bin/bash
        mkdir -p $outputdir
        echo -e "\n***** SPLIT VCFs for POP: ${pop}, VCF file: $vcf_file"
        #bash "${launchDir}/bin/step0_split_vcfs_parallel.sh" \
        #    --vcf_path "$vcf_file" \
        #    --output_dir "$outputdir" \
        #    --num_of_cores 20
        bash "${launchDir}/bin/step0_split_vcfs_parallel.sh" \
            -v "$vcf_file" \
            -d "$outputdir" \
            -c 20
        """
}


process filter_vcf_generegions_withgnomad {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file), val(gene_regions_file), val(gnomad_raw_vcf_file) 

    output:
        path "region_based_vcfs_withgnomad/${pop}/${pop}_filtered.vcf.gz", emit: region_filtered_vcf
        path "region_based_vcfs_withgnomad/${pop}/gnomad_filtered.vcf.gz", emit: region_filtered_gnomadvcf
        path "region_based_vcfs_withgnomad/${pop}/${pop}_filtered.vcf.gz.tbi", emit: filt_vcf_index
        path "region_based_vcfs_withgnomad/${pop}/gnomad_filtered.vcf.gz.tbi", emit: filt_gnomadvcf_index

    script:
        def outputdir = "region_based_vcfs_withgnomad/${pop}"

        """
        #!/usr/bin/bash
        mkdir -p $outputdir
        echo -e "\n***** FILTER RAREVARIANT for POP: ${pop}, VCF file: $vcf_file, GENE region file: $gene_regions_file, GNOMAD VCF file: $gnomad_raw_vcf_file"
        bash "${launchDir}/bin/step0.1_filter_vcf_generegions.sh" \
            -o "$outputdir" \
            -v "$vcf_file" \
            -r "$gene_regions_file" \
            -g "$gnomad_raw_vcf_file" \
            -c "$pop"
        """
}


process filter_vcf_generegions {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file), val(gene_regions_file) 

    output:
        path "region_based_vcfs/${pop}/*_filtered.vcf.gz", emit: region_filtered_vcf
        path "region_based_vcfs/${pop}/*_filtered.vcf.gz.tbi", emit: filt_vcf_index

    script:
        def outputdir = "region_based_vcfs/${pop}"

        """
        #!/usr/bin/bash
        mkdir -p $outputdir
        echo -e "\n***** FILTER RAREVARIANT for POP: ${pop}, VCF file: $vcf_file, GENE region file: $gene_regions_file"
        bash "${launchDir}/bin/step0.1_filter_vcf_generegions-vcfonly.sh" \
            -o "$outputdir" \
            -v "$vcf_file" \
            -r "$gene_regions_file" \
            -c "$pop"
        """
}


process subset_pop_and_qc {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file), val(pop_list) 

    output:
        path "subset_pop_vcf/${pop}/*_subset.vcf.gz", emit: subset_vcf_ch
        path "subset_pop_vcf/${pop}/*_subset.vcf.gz.tbi", emit: subset_vcf_index_ch

    script:
        def outputdir = "subset_pop_vcf/${pop}"

        """
        #!/usr/bin/bash
        mkdir -p $outputdir
        echo -e "\n***** FILTER POP SUBSET: ${pop}, VCF file: $vcf_file, POP Subset List file: $pop_list"
        bash "${launchDir}/bin//step0.2_subset_pop_and_qc.sh" \
            -o "$outputdir" \
            -v "$vcf_file" \
            -c "$pop" \
            -p "$pop_list"
        """
}


process filter_rare_variants_by_maf {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(subset_vcf_file) 

    output:
        path "rare_vars/${pop}/*_rare.vcf.gz", emit: rare_variants_vcf_ch
        path "rare_vars/${pop}/*_rare.vcf.gz.tbi", emit: rare_variants_index_ch
        path "rare_vars/${pop}/*_rare.bed", emit: rare_variants_bed_ch

    script:
        def outputdir = "rare_vars/${pop}"

        """
        #!/usr/bin/bash
        mkdir -p $outputdir
        echo -e "\n***** RARE VARIANTS CALL: ${pop}, VCF file: $subset_vcf_file"
        bash "${launchDir}/bin/step0.3_filter_rare_variants_by_maf.sh" \
            -o "$outputdir" \
            -v "$subset_vcf_file" \
            -c "$pop"
        """
}


process concat_rare_variants_vcf_awk {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        val(pop)
        val(rarevar_chr_vcflist) 

    output:
        path "rare_vars/${pop}/concatenated/*_rare.vcf.gz", emit: concat_rare_variants_vcf_ch

    script:
        def outputdir = "rare_vars/${pop}/concatenated"
        def rarevar_chr_vcflist_string = rarevar_chr_vcflist.join(' ')
        def cohort_rare_vcf="${outputdir}/${rarevar_chr_vcflist[1].Name.replaceFirst(/chr[0-9XY]+_/, '')}" //replaceFirst also avail

        """
        num_threads=\$(nproc)
        mkdir -p $outputdir
        echo -e "\n***** CONCATENATED RARE VARIANTS VCF POP: ${pop}, Number of threads: \$num_threads, RV VCF FILE LIST: $rarevar_chr_vcflist, OUTPUT FILE NAME 1: $cohort_rare_vcf"
        echo "${rarevar_chr_vcflist_string}" | tr ' ' '\n' | sort -V > ${outputdir}/unsorted_file_list.txt

        # Sort the file list based on the chromosome number extracted from the filename
         awk -F'_' '
         {
              match(\$0, /chr([0-9XY]+)_filtered_subset_rare/, arr);
              chr = arr[1];
              if (chr == "X") chr = 23;
              else if (chr == "Y") chr = 24;
              print chr, \$0;

        }' ${outputdir}/unsorted_file_list.txt | sort -k1,1n | cut -d' ' -f2- > ${outputdir}/sorted_file_list.txt

        # Concatenate VCFs based on file lists
        bcftools concat --file-list ${outputdir}/sorted_file_list.txt --output $cohort_rare_vcf --output-type z --threads \$num_threads

        echo "**** Indexing VCF"
        bcftools index --tbi $cohort_rare_vcf --threads \$num_threads
        """
}


process concat_rare_variants_vcf {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        val(pop)
        val(rarevar_chr_vcflist) 

    output:
        path "rare_vars/${pop}/concatenated/*_rare.vcf.gz", emit: concat_rare_variants_vcf_ch

    script:
        def outputdir = "rare_vars/${pop}/concatenated"
        def rarevar_chr_vcflist_string = rarevar_chr_vcflist.join(' ')
        def cohort_rare_vcf="${outputdir}/${rarevar_chr_vcflist[1].Name.replaceFirst(/chr[0-9XY]+_/, '')}" //replaceFirst also avail

        """
        # Set number of threads
        #num_threads=\$((\$(nproc) / 4))
        num_threads=\$(nproc)
        echo "Using \$num_threads threads for parallel processing."
         
        mkdir -p $outputdir
        echo -e "\n***** CONCATENATED RARE VARIANTS VCF POP: ${pop}, Number of threads: \$num_threads, RV VCF FILE LIST: $rarevar_chr_vcflist, OUTPUT FILE NAME 1: $cohort_rare_vcf"
        echo "${rarevar_chr_vcflist_string}" | tr ' ' '\n' | sort -V > ${outputdir}/unsorted_file_list.txt

        # Extract basenames, sort them, and map back to full paths
        awk '{print \$0}' ${outputdir}/unsorted_file_list.txt | while read line; do echo "\$line \$(basename \$line)"; done | sort -k2,2V | awk '{print \$1}' > ${outputdir}/sorted_final_file_list.txt

        # Concatenate VCFs based on file lists
        bcftools concat --file-list ${outputdir}/sorted_final_file_list.txt --output $cohort_rare_vcf --output-type z --threads \$num_threads

        echo "**** Indexing VCF"
        bcftools index --tbi $cohort_rare_vcf --threads \$num_threads
        """
}


process concat_rare_variants_bed {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        val(pop)
        val(rarevar_chr_bedlist) 

    output:
        //path "rare_vars/${pop}/concatenated/*_rare.bed", emit: rare_variants_bed_ch
        path "rare_vars/${pop}/concatenated/*_rare.bed", emit: concat_rare_variants_bed_ch

    script:
        def outputdir = "rare_vars/${pop}/concatenated"
        def rarevar_chr_bedlist_string = rarevar_chr_bedlist.join(' ')
        def cohort_rare_bed="${outputdir}/${rarevar_chr_bedlist[1].Name.replaceFirst(/chr[0-9XY]+_/, '')}" //replaceFirst also avail

        //def cohort_rare_bed2="${outputdir}/${rarevar_chr_bedlist[1].baseName}" //.bed removed
        //def cohort_rare_bed3="${outputdir}/${rarevar_chr_bedlist[1].simpleName}" //all multiple "." suffix removed

        """
        mkdir -p $outputdir
        echo -e "\n***** CONCATENATED RARE VARIANTS BED: ${pop}, RV BED FILE LIST: $rarevar_chr_bedlist, RV BED STRING LIST: $rarevar_chr_bedlist_string, OUTPUT FILE NAME 1: $cohort_rare_bed"
        echo "${rarevar_chr_bedlist_string}" > ${outputdir}/file_list.txt
        cat ${outputdir}/file_list.txt | xargs cat | bedtools sort -i - > "${cohort_rare_bed}"
        """
}


process identify_gnomad_common_variants {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(gnomad_chr_vcf), val(cohort_rare_bed), val(ancestry)

    output:
        path "gnomad_common_vars/${pop}/*_common.vcf.gz", emit: common_variants_ch
        path "gnomad_common_vars/${pop}/*_common.vcf.gz.tbi", emit: common_variants_index_ch
        path "gnomad_common_vars/${pop}/*_common.bed", emit: common_variants_bed_ch

    script:
        def outputdir = "gnomad_common_vars/${pop}"

        """
        #echo_test
        mkdir -p $outputdir
        echo -e "\n***** RARE VARIANTS CALL: ${pop}, GNOMAD CHROM VCF file: $gnomad_chr_vcf, RARE Bed file: $cohort_rare_bed, Ancestry: $ancestry"
        bash "${launchDir}/bin/step0.4_identify_gnomad_common_variants.sh" \
            -o "$outputdir" \
            -g "$gnomad_chr_vcf" \
            -r "$cohort_rare_bed" \
            -c "$ancestry"
        """
}


process identify_gnomad_common_variants_chromwise {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(chrom), val(gnomad_chr_vcf), val(cohort_rare_bed), val(pop), val(ancestry)

    output:
        path "gnomad_common_vars_chromwise/${pop}/*_common.vcf.gz", emit: common_variants_ch
        path "gnomad_common_vars_chromwise/${pop}/*_common.vcf.gz.tbi", emit: common_variants_index_ch
        path "gnomad_common_vars_chromwise/${pop}/*_common.bed", emit: common_variants_bed_ch

    script:
        def outputdir = "gnomad_common_vars_chromwise/${pop}"

        """
        #echo_test
        mkdir -p $outputdir
        echo -e "\n***** RARE VARIANTS CALL for chrom: $chrom, POP:  ${pop}, GNOMAD CHROM VCF file: $gnomad_chr_vcf, RARE Bed file: $cohort_rare_bed, Ancestry: $ancestry"
        bash "${launchDir}/bin/step0.4_identify_gnomad_common_variants.sh" \
            -o "$outputdir" \
            -g "$gnomad_chr_vcf" \
            -r "$cohort_rare_bed" \
            -c "$ancestry"
        """
}


process concat_gnomad_common_varaints_bed {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        val(pop)
        val(gnomadcommon_chr_bedlist) 

    output:
        path "gnomad_common_vars_chromwise/${pop}/concatenated/*_common.bed", emit: concat_gnomadcommon_variants_bed_ch

    script:
        def outputdir = "gnomad_common_vars_chromwise/${pop}/concatenated"
        def gnomadcommon_chr_bedlist_string = gnomadcommon_chr_bedlist.join(' ')
        def cohort_gnomadcommon_bed="${outputdir}/${gnomadcommon_chr_bedlist[1].Name.replaceFirst(/chr[0-9XY]+_/, '')}" //replaceFirst also avail

        """
        mkdir -p $outputdir
        echo -e "\n***** CONCATENATED GNOMAD COMMON VARIANTS BED: ${pop}, RV BED FILE LIST: $gnomadcommon_chr_bedlist, RV BED STRING LIST: $gnomadcommon_chr_bedlist_string, OUTPUT FILE NAME 1: $cohort_gnomadcommon_bed"
        echo "${gnomadcommon_chr_bedlist_string}" > ${outputdir}/file_list.txt
        cat ${outputdir}/file_list.txt | xargs cat | bedtools sort -i - > "${cohort_gnomadcommon_bed}"
        """
}


process rarevars_map_to_genes {
    debug true
    tag { "JOB ${pop}" }
    //label 'watershed_rv_env'
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(gnomad_common_bed), val(cohort_rare_vcf), val(gene_regions_bed)

    output:
        path "rv_gene_map/${pop}/cohort_${pop}_rareQC.vcf.gz", emit: cohort_rare_qc_vcf
        path "rv_gene_map/${pop}/cohort_${pop}_rareQC.indiv.txt", emit: indiv_at_rv
        path "rv_gene_map/${pop}/gene-${pop}-rv.raw.txt", emit: rv_sites_raw

    script:
        def outputdir = "rv_gene_map/${pop}"

        """
        mkdir -p $outputdir
        echo -e "\n***** RARE VARIANTS CALL for POP:  ${pop}, GNOMAD CHROM BED file: $gnomad_common_bed, RARE VCF file: $cohort_rare_vcf"
        bash "${launchDir}/bin/step0.5_rarevars_map_to_genes.sh" \
            -o "$outputdir" \
            -g "${gnomad_common_bed}" \
            -r "${cohort_rare_vcf}" \
            -b "${gene_regions_bed}" \
            -c "${pop}"
        """
}


process rarevars_map_to_gene_indivs_format {
    //build_gene_indiv_rare_variants
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/0_rv_files//", mode: 'copy'

    input:
        tuple val(pop), path(rv_sites_raw)

    output:
        path "rv_gene_map/${pop}/gene-${pop}-rv.txt", emit: gene_indiv_rare_variants

    script:
        def outputdir = "rv_gene_map/${pop}"
        def outputfile = "$outputdir/gene-${pop}-rv.txt"

        """
        mkdir -p $outputdir
        echo -e "\n**** Rare variants per each gene-individual pair for POP: $pop, RV SITES RAW file: $rv_sites_raw"
        Rscript "${launchDir}/bin/step0.6_gene_indiv_rare_variants.R" \
        --rv_sites "${rv_sites_raw}" \
        --popname "${pop}" \
        --outfile "${outputfile}"

        echo "Rare variants gene-individual pair file generated. Check result file: ${outputfile}"
        """
}


process convert_rare_variants_to_bed {
    debug true
    tag { "JOB ${pop}" }
    label 'watershed_r_env'
    publishDir "$params.out_dir/outputs/0_rv_files//", mode: 'copy'

    input:
        tuple val(pop), path(rv_sites)

    output:
        path "rv_gene_map/${pop}/gene-${pop}-rv.bed", emit: gene_indiv_rare_variants_bed

    script:
        def outputdir = "rv_gene_map/${pop}"
        def outputfile = "$outputdir/gene-${pop}-rv.bed"

        """
        mkdir -p $outputdir
        echo -e "\n**** Rare variants per each gene-individual pair for POP: $pop, RV SITES file: $rv_sites"
        Rscript "${launchDir}/bin/step8.1_rare_variants_to_bed.R" \
        --rv_sites "${rv_sites}" \
        --outfile "${outputfile}"

        echo "Bed conversion of rare varaints. Check result file: ${outputfile}"
        echo -e "\n\n******Rare Variant Operation DONE!*******\n"
        """
}


workflow CALL_RAREVARIANTS {
    take:
    pops_ch
    rv_file_ch
    gencode_region_file_ch
    gnomad_file_ch
    subjids_pops_list_ch
    gnomad_generegion_filelist_ch
    ancestry_ch

    main:
    if (params.analysis == "GLOBAL") {
        
        //Channel.value("$params.cohort_name").set { pops_ch }
        //pops_ch.merge(rv_file_ch, vep_loftee_file_ch)

        // Combine the final outputs into a single channel
        combined_splitvcf_ch = pops_ch
            .merge(rv_file_ch)
            .map{ pops, rv_file  -> tuple(pops, rv_file)
        }

        combined_splitvcf_ch.view({ "SUBWORKFLOW RV FILE combined output: $it" })

        // Split VCFs by chromosomes
        splitvcfs_by_chrom(combined_splitvcf_ch)

        // Split files output
        vcffiles_ch = splitvcfs_by_chrom.out.splitvcfs

        vcffiles_ch.flatten()
            .take(24)
            .set { vcf_channel }


        /*
        vcffiles_ch.flatten()
            .filter { vcf -> !vcf.name.contains("chr1") && !vcf.name.contains("chr2.") }
            .take(3)
            .set { vcf_channel }

        // Combine the final outputs into a single channel
        pops_ch
            .combine(vcf_channel)
            .combine(gencode_region_file_ch)
            .combine(gnomad_file_ch)
            .set { combined_gnomad_ch }

        combined_gnomad_ch.view({"Here's the RARE VARIANT VCF COMBINED INPUT CHANNEL: $it"})

        // Filter gene regions from the vcfs plus gnomads
        filter_vcf_generegions_withgnomad(combined_gnomad_ch)
        */



        // Combine the final outputs into a single channel
        pops_ch
            .combine(vcf_channel)
            .combine(gencode_region_file_ch)
            .set { combined_vcf_ch }

        combined_vcf_ch.view({"Here's the RARE VARIANT VCF COMBINED INPUT CHANNEL: $it"})

        // Filter gene regions from the vcfs with parallel process
        filter_vcf_generegions(combined_vcf_ch)


        
        // Capture gene region vcf files from 22 chroms parallel process
        vcf_generegion_file_ch = filter_vcf_generegions.out.region_filtered_vcf

        // Combine necessary channels to prepare for subsetting population VCF
        pops_ch
            .combine(vcf_generegion_file_ch)
            .combine(subjids_pops_list_ch)
            .set { combined_popsubset_ch }

        combined_popsubset_ch.view({"Here's the INDIV POP SUBSET VCF COMBINED INPUT CHANNEL: $it"})

        // Subset population vcf 
        subset_pop_and_qc(combined_popsubset_ch)

        // Capture pop subset vcf files from 22 chroms parallel process
        pop_subset_vcf_ch = subset_pop_and_qc.out.subset_vcf_ch



        // Combine necessary channels to prepare for subsetting population VCF
        pops_ch
            .combine(pop_subset_vcf_ch)
            .set { combined_rarevar_ch }

        combined_rarevar_ch.view({"Here's the INDIV POP SUBSET VCF COMBINED INPUT CHANNEL: $it"})

        // Filter rare variants 
        filter_rare_variants_by_maf(combined_rarevar_ch)

        // Rare bed files
        rare_bedfilelist_ch = filter_rare_variants_by_maf.out.rare_variants_bed_ch.collect()
        rare_vcffilelist_ch = filter_rare_variants_by_maf.out.rare_variants_vcf_ch.collect()


        //rare_bedfilelist_ch.view({ "Here's the RARE BED REGION filelist: $it" })
        //gnomad_generegion_filelist_ch.view({ "Here's the gnomAD region filelist: $it" })

       
        concat_rare_variants_bed(pops_ch, rare_bedfilelist_ch)
        rarebed_generegion_file_ch = concat_rare_variants_bed.out.concat_rare_variants_bed_ch

        pops_ch
            .combine(gnomad_generegion_filelist_ch)
            .combine(rarebed_generegion_file_ch)
            .combine(ancestry_ch)
            .set { combined_gnomadcommon_ch }


        combined_gnomadcommon_ch.view({ "Here's the COMBINED COMMON GNOMAD filelist: $it" })
        /*identify_gnomad_common_variants(combined_gnomadcommon_ch)*/

        /*gnomad_common_ch = identify_gnomad_common_variants.out.common_variants_bed_ch.collect()*/
        /*gnomad_common_ch.view({ "Here's the COLLECT VERSION OF COMMON GNOMAD BED filelist: $it" })*/
        
        // Function to extract chromosome number from file paths
        def extractChroms = { filePath ->
            def match = filePath =~ /(chr(\d+|X|Y))/
            return match ? match[0][1] : null
        }

        // For capturing only num: "1", "2", "X" instead of "chr1", "chr2", "chrX" above     
        /* /chr(\d+|X|Y)/ */

        // Create a key-value pair channel for rare_variant_ch
        rare_bedfilelist_ch.flatten().map { file ->
            def chrom = extractChroms(file.toString())
            return [chrom, file]
        }.set { rare_variant_keyed_ch }

        // Create a key-value pair channel for gnomad_generegion_filelist_ch
        gnomad_generegion_filelist_ch.map { file ->
            def chrom = extractChroms(file.toString())
            return [chrom, file]
        }.set { gnomad_generegion_keyed_ch }


        rare_variant_keyed_ch.view({ "Here's the OUTPUT OF rare_variant_keyed_ch: $it" })
        gnomad_generegion_keyed_ch.view({ "Here's the OUTPUT OF gnomad_generegion_keyed_ch: $it" }) 
        
        // Joining channels based on the key
        gnomad_generegion_keyed_ch
            .join(rare_variant_keyed_ch)
            .merge(pops_ch)
            .merge(ancestry_ch)
            .map { chrom, gnomad_filtered, rare_variant, pops_ch, ancestry_ch ->
                tuple(chrom, gnomad_filtered, rare_variant, pops_ch, ancestry_ch)
            }.set { combined_keyed_ch }

        combined_keyed_ch.view({ "Here's the OUTPUT OF JOINED KEY CHANNEL: $it" }) 
        identify_gnomad_common_variants_chromwise(combined_keyed_ch)

        // Parameters
        identify_gnomad_common_variants_chromwise
            .out
            .common_variants_bed_ch
            .set{ gnomad_common_chromwise_ch }


        /*    .collect() */

        //Channel.fromPath("$params.gnomad_common_chromwise_filelist").set { gnomad_common_chromwise_ch }
        gnomad_common_chromwise_ch.view({ "Here's the OUTPUT OF GNOMAD CHROMWISE OPERATION: $it" })

        // Concatenate gnomad common bed files     
        concat_gnomad_common_varaints_bed(pops_ch, gnomad_common_chromwise_ch.collect())

        concat_gnomadcommon_bed_file_ch = concat_gnomad_common_varaints_bed.out.concat_gnomadcommon_variants_bed_ch
        gnomad_common_chromwise_ch.view({ "Here's the CONCAT OUTPUT OF GNOMAD COMMON BED FILE: $it" })


        // Concatenate Rare variant vcf files        
        concat_rare_variants_vcf(pops_ch, rare_vcffilelist_ch)
        concat_rare_variants_vcf
            .out
            .concat_rare_variants_vcf_ch
            .set { concat_rare_variants_vcf_ch }

        concat_rare_variants_vcf_ch.view({ "Here's the CONCATENATED RARE VARIANT VCF FILE: $it" })

        // Map rare variants to genes
        pops_ch
            .merge(concat_gnomadcommon_bed_file_ch)
            .merge(concat_rare_variants_vcf_ch)
            .merge(gencode_region_file_ch)
            .set { combined_rarevars_genemap_ch }

        combined_rarevars_genemap_ch.view({ "Here's the COMBINED INPUT TUPLE for RARE VARIANT TO GENE MAP: $it" })
        
        //Run
        rarevars_map_to_genes(combined_rarevars_genemap_ch)
        rarevars_map_to_genes
            .out
            .rv_sites_raw
            .set {gene_pop_raw_rv_file_ch }

        // Map rare variants to gene-indivs formatted(rearranged)
        pops_ch
            .merge(gene_pop_raw_rv_file_ch)
            .set { combined_rarevars_geneindiv_map_ch }

        combined_rarevars_geneindiv_map_ch.view({ "Here's the COMBINED INPUT TUPLE for RARE VARIANT TO GENE INDIV MAP REFORMATTING: $it" })
        
        //Run
        rarevars_map_to_gene_indivs_format(combined_rarevars_geneindiv_map_ch)
        rarevars_map_to_gene_indivs_format
            .out
            .gene_indiv_rare_variants
            .set { gene_indiv_rv_ch }

        // Convert rare variants to bed
        pops_ch
            .merge(gene_indiv_rv_ch)
            .set { combined_rarevars_geneindiv_bed_ch }

        combined_rarevars_geneindiv_bed_ch.view({ "Here's the COMBINED INPUT TUPLE for RARE VARIANT TO BED CONVERSION: $it" })
        
        //Run
        convert_rare_variants_to_bed(combined_rarevars_geneindiv_bed_ch)
        convert_rare_variants_to_bed
            .out
            .gene_indiv_rare_variants_bed
            .set { gene_indiv_rv_bed_ch }



        //.view({ chrom, gnomad_filtered, rare_variant -> "Final combined tuple: ($chrom, $gnomad_filtered, $rare_variant)" })
        //.filter { key, lognorm_expr, peer_factors, covariates -> key == "EUR" || key == "AFR" }
        /*
        pops_ch
            .merge(vcf_generegion_file_ch, subjids_pops_list_ch)
            .map { pops, vcf_generegion_file, pops_list -> tuple(pops, vcf_generegion_file, pops_list) }
            .set { combined_popsubset_ch }
       */


       /*------------------------------------------------------------------*/

        /*
        vcffiles_ch.flatten()
            .filter { vcf -> !vcf.name.contains("chr1") && !vcf.name.contains("chr2.") }
            .take(3)
            .combine(pops_ch)
            .combine(gencode_region_file_ch)
            .combine(gnomad_file_ch)
            .set { vcf_comb_channel }

        vcf_comb_channel.view({"Here's the RARE VARIANT VCF COMBINED INPUT CHANNEL: $it"})
        */


         //vcffiles_ch.view({"Here's the VCF file example: $it[0]"})
        //.take(params.num_chromosomes)
        //Run Rare variant process as a separate process
        //call_cadd_annotation(vcf_channel)

       /* vcffiles_ch.flatten()
            .filter { vcf -> !vcf.name.contains("chr1") && !vcf.name.contains("chr2.") }
            .take(3)
            .map { vcf -> tuple("$params.cohort_name", vcf) }
            .set { vcf_channel }
        */


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
    splitvcf_files = splitvcfs_by_chrom.out.splitvcfs
    splitvcf_files_flat = splitvcfs_by_chrom.out.splitvcfs.flatten().take(4)
    filt_vcf_file = filter_vcf_generegions.out.region_filtered_vcf
    //filt_gnomadvcf_file = filter_vcf_generegions.out.region_filtered_gnomadvcf
    pop_subset_vcf_file = subset_pop_and_qc.out.subset_vcf_ch
    rarevar_vcf_file = filter_rare_variants_by_maf.out.rare_variants_vcf_ch
    rarevar_bed_file = filter_rare_variants_by_maf.out.rare_variants_bed_ch
    /*gnomad_common_bed_file = identify_gnomad_common_variants.out.common_variants_bed_ch.collect()*/
    gnomad_common_bed_chrom_file = identify_gnomad_common_variants_chromwise.out.common_variants_bed_ch.flatten().take(4)
    concat_rare_variants_vcf_file = concat_rare_variants_vcf.out.concat_rare_variants_vcf_ch
    gene_pop_raw_rv_file = rarevars_map_to_genes.out.rv_sites_raw
    gene_indiv_rv_file = gene_indiv_rv_ch
    gene_indiv_rv_bed_file = gene_indiv_rv_bed_ch

    /*splitvcf_files_flat = splitvcfs_by_chrom.out.splitvcfs.flatten()
                    .filter { vcf -> !vcf.name.contains("chr1") && !vcf.name.contains("chr2.") }
                    .take(3)
    */
}

