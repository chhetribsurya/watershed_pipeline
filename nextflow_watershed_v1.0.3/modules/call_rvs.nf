#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process splitvcfs_by_chrom {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file)

    output:
        path "splitvcfs/${pop}/*.vcf.gz", emit: splitvcfs
        path "splitvcfs/${pop}/*.vcf.gz.tbi", emit: indexes

    script:
        def outputdir = "splitvcfs/${pop}"

        """
        mkdir -p $outputdir
        bash "${launchDir}/bin/step0_split_vcfs_parallel.sh" \
            -v "$vcf_file" \
            -d "$outputdir" \
            -c 20
        """
}


process filter_vcf_generegions_withgnomad {
    tag { "JOB ${pop}" }
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
        mkdir -p $outputdir
        bash "${launchDir}/bin/step0.1_filter_vcf_generegions.sh" \
            -o "$outputdir" \
            -v "$vcf_file" \
            -r "$gene_regions_file" \
            -g "$gnomad_raw_vcf_file" \
            -c "$pop"
        """
}


process filter_vcf_generegions {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file), val(gene_regions_file) 

    output:
        path "region_based_vcfs/${pop}/*_filtered.vcf.gz", emit: region_filtered_vcf
        path "region_based_vcfs/${pop}/*_filtered.vcf.gz.tbi", emit: filt_vcf_index

    script:
        def outputdir = "region_based_vcfs/${pop}"

        """
        mkdir -p $outputdir
        bash "${launchDir}/bin/step0.1_filter_vcf_generegions-vcfonly.sh" \
            -o "$outputdir" \
            -v "$vcf_file" \
            -r "$gene_regions_file" \
            -c "$pop"
        """
}


process subset_pop_and_qc {
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        tuple val(pop), val(vcf_file), val(pop_list) 

    output:
        path "subset_pop_vcf/${pop}/*_subset.vcf.gz", emit: subset_vcf_ch
        path "subset_pop_vcf/${pop}/*_subset.vcf.gz.tbi", emit: subset_vcf_index_ch

    script:
        def outputdir = "subset_pop_vcf/${pop}"

        """
        mkdir -p $outputdir
        bash "${launchDir}/bin//step0.2_subset_pop_and_qc.sh" \
            -o "$outputdir" \
            -v "$vcf_file" \
            -c "$pop" \
            -p "$pop_list"
        """
}


process filter_rare_variants_by_maf {
    tag { "JOB ${pop}" }
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
        mkdir -p $outputdir
        bash "${launchDir}/bin/step0.3_filter_rare_variants_by_maf.sh" \
            -o "$outputdir" \
            -v "$subset_vcf_file" \
            -c "$pop"
        """
}


process concat_rare_variants_vcf_awk {
    tag { "JOB ${pop}" }
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
        bcftools index --tbi $cohort_rare_vcf --threads \$num_threads
        """
}


process concat_rare_variants_vcf {
    tag { "JOB ${pop}" }
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
    tag { "JOB ${pop}" }
    publishDir "$params.out_dir/outputs/0_rv_files/", mode : 'copy'

    input:
        val(pop)
        val(rarevar_chr_bedlist) 

    output:
        path "rare_vars/${pop}/concatenated/*_rare.bed", emit: concat_rare_variants_bed_ch

    script:
        def outputdir = "rare_vars/${pop}/concatenated"
        def rarevar_chr_bedlist_string = rarevar_chr_bedlist.join(' ')
        def cohort_rare_bed="${outputdir}/${rarevar_chr_bedlist[1].Name.replaceFirst(/chr[0-9XY]+_/, '')}" //replaceFirst also avail


        """
        mkdir -p $outputdir
        echo "${rarevar_chr_bedlist_string}" > ${outputdir}/file_list.txt
        cat ${outputdir}/file_list.txt | xargs cat | bedtools sort -i - > "${cohort_rare_bed}"
        """
}


process identify_gnomad_common_variants {
    tag { "JOB ${pop}" }
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
        mkdir -p $outputdir
        bash "${launchDir}/bin/step0.4_identify_gnomad_common_variants.sh" \
            -o "$outputdir" \
            -g "$gnomad_chr_vcf" \
            -r "$cohort_rare_bed" \
            -c "$ancestry"
        """
}


process identify_gnomad_common_variants_chromwise {
    tag { "JOB ${pop}" }
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
        mkdir -p $outputdir
        bash "${launchDir}/bin/step0.4_identify_gnomad_common_variants.sh" \
            -o "$outputdir" \
            -g "$gnomad_chr_vcf" \
            -r "$cohort_rare_bed" \
            -c "$ancestry"
        """
}


process concat_gnomad_common_varaints_bed {
    tag { "JOB ${pop}" }
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
        echo "${gnomadcommon_chr_bedlist_string}" > ${outputdir}/file_list.txt
        cat ${outputdir}/file_list.txt | xargs cat | bedtools sort -i - > "${cohort_gnomadcommon_bed}"
        """
}


process rarevars_map_to_genes {
    tag { "JOB ${pop}" }
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
        bash "${launchDir}/bin/step0.5_rarevars_map_to_genes.sh" \
            -o "$outputdir" \
            -g "${gnomad_common_bed}" \
            -r "${cohort_rare_vcf}" \
            -b "${gene_regions_bed}" \
            -c "${pop}"
        """
}


process rarevars_map_to_gene_indivs_format {
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
        Rscript "${launchDir}/bin/step0.6_gene_indiv_rare_variants.R" \
        --rv_sites "${rv_sites_raw}" \
        --popname "${pop}" \
        --outfile "${outputfile}"
        """
}


process convert_rare_variants_to_bed {
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
        Rscript "${launchDir}/bin/step8.1_rare_variants_to_bed.R" \
        --rv_sites "${rv_sites}" \
        --outfile "${outputfile}"
        """
}

