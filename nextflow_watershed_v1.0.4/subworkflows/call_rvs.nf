#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
  splitvcfs_by_chrom;
  filter_vcf_generegions;
  subset_pop_and_qc;
  filter_rare_variants_by_maf;
  concat_rare_variants_bed;
  identify_gnomad_common_variants_chromwise;
  concat_gnomad_common_varaints_bed;
  concat_rare_variants_vcf;
  rarevars_map_to_genes;
  rarevars_map_to_gene_indivs_format;
  convert_rare_variants_to_bed
} from "../modules/call_rvs" params(params)

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

        // Combine the final outputs into a single channel
        combined_splitvcf_ch = pops_ch
            .merge(rv_file_ch)
            .map{ pops, rv_file  -> tuple(pops, rv_file)
        }

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

        // Filter gene regions from the vcfs with parallel process
        filter_vcf_generegions(combined_vcf_ch)

        // Capture gene region vcf files from 22 chroms parallel process
        vcf_generegion_file_ch = filter_vcf_generegions.out.region_filtered_vcf

        // Combine necessary channels to prepare for subsetting population VCF
        pops_ch
            .combine(vcf_generegion_file_ch)
            .combine(subjids_pops_list_ch)
            .set { combined_popsubset_ch }

        // Subset population vcf 
        subset_pop_and_qc(combined_popsubset_ch)

        // Capture pop subset vcf files from 22 chroms parallel process
        pop_subset_vcf_ch = subset_pop_and_qc.out.subset_vcf_ch

        // Combine necessary channels to prepare for subsetting population VCF
        pops_ch
            .combine(pop_subset_vcf_ch)
            .set { combined_rarevar_ch }

        // Filter rare variants 
        filter_rare_variants_by_maf(combined_rarevar_ch)

        // Rare bed files
        rare_bedfilelist_ch = filter_rare_variants_by_maf.out.rare_variants_bed_ch.collect()
        rare_vcffilelist_ch = filter_rare_variants_by_maf.out.rare_variants_vcf_ch.collect()

        concat_rare_variants_bed(pops_ch, rare_bedfilelist_ch)
        rarebed_generegion_file_ch = concat_rare_variants_bed.out.concat_rare_variants_bed_ch

        pops_ch
            .combine(gnomad_generegion_filelist_ch)
            .combine(rarebed_generegion_file_ch)
            .combine(ancestry_ch)
            .set { combined_gnomadcommon_ch }


        // Function to extract chromosome number from file paths
        def extractChroms = { filePath ->
            def match = filePath =~ /(chr(\d+|X|Y))/
            return match ? match[0][1] : null
        }

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

        
        // Joining channels based on the key
        gnomad_generegion_keyed_ch
            .join(rare_variant_keyed_ch)
            .merge(pops_ch)
            .merge(ancestry_ch)
            .map { chrom, gnomad_filtered, rare_variant, pops_ch, ancestry_ch ->
                tuple(chrom, gnomad_filtered, rare_variant, pops_ch, ancestry_ch)
            }.set { combined_keyed_ch }

        identify_gnomad_common_variants_chromwise(combined_keyed_ch)

        // Parameters
        identify_gnomad_common_variants_chromwise
            .out
            .common_variants_bed_ch
            .set{ gnomad_common_chromwise_ch }


        // Concatenate gnomad common bed files     
        concat_gnomad_common_varaints_bed(pops_ch, gnomad_common_chromwise_ch.collect())

        concat_gnomadcommon_bed_file_ch = concat_gnomad_common_varaints_bed.out.concat_gnomadcommon_variants_bed_ch


        // Concatenate Rare variant vcf files        
        concat_rare_variants_vcf(pops_ch, rare_vcffilelist_ch)
        concat_rare_variants_vcf
            .out
            .concat_rare_variants_vcf_ch
            .set { concat_rare_variants_vcf_ch }

        // Map rare variants to genes
        pops_ch
            .merge(concat_gnomadcommon_bed_file_ch)
            .merge(concat_rare_variants_vcf_ch)
            .merge(gencode_region_file_ch)
            .set { combined_rarevars_genemap_ch }

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

        //Run
        convert_rare_variants_to_bed(combined_rarevars_geneindiv_bed_ch)
        convert_rare_variants_to_bed
            .out
            .gene_indiv_rare_variants_bed
            .set { gene_indiv_rv_bed_ch }


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

