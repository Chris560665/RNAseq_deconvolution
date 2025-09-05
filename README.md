# RNAseq_deconvolution
1. The preprocess of a gzvcf file from a mixture transcriptome sequencing data could be conducted by
            ./mixture_preprocess.sh -s sample_basename
   in which "All.chr.EAS.recode.vcf.gz" is accessed from the 1000G project and obtained by merging vcf files from all chromosomes and extracting peoples from East Asian populations using vcftools, RNAedit_pos.txt contained the positions of SNPs with editing events accessed from REDIportal, and "imprintgene_pos.txt" contained the position of imprint gene.
2. The SNP deconvolution for two-body-fluid mixture should be performed through the following code:
            ./mixture_phenotype.sh -s sample_basename -cF contributed_body_fluid_1 -cS contributed_body_fluid_2 -bF bed_file_1 -bS bed_file_2
            ./phenotype_frequency.sh -s sample_basename -cF contributed_body_fluid_1 -cS contributed_body_fluid_2 -bF bed_file_1 -bS bed_file_2
4. The bed files used in the scripts should include the information of chromosome, start position, end position and gene symbol of relative specific genes, for example, for VB-SE mixtures, "bed_file_1" should including the above information of overexpressed genes in venous blood compared to semen, and bed_file_2 should including the above information of overexpressed genes in semen compared to blood. 
7. The python scripts "summarized_data.py" were used to automatically assign the mixture phenotype to a phenotype frequency, which could be implemented by the code (a example for the results of VB contributor):
   python summarized_data.py -m sample_basename.contributed_body_fluid_1.mix_pheno.txt -e sample_basename.contributed_body_fluid_1_EAS_phenofreq.txt -o sample_basename.contributed_body_fluid_1.xlsx
in this code, the "sample_basename.contributed_body_fluid_1.mix_pheno.txt" is the result of VB contributor accessed by "mixture_phenotype.sh", and the "sample_basename.contributed_body_fluid_1_EAS_phenofreq.txt" is the result of VB contributor accessed by "phenotype_frequency.sh".
