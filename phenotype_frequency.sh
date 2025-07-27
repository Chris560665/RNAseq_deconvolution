#This code is exemplified for the deconvolution of VB-SE mixtures, but could also be extended to other multi-body-fluid-mixtures.
#The folders' name could be revised, and the folders(BF_1000G) and some of corresponding codes could be added according to the deconvoluted mixture types.
#The bed files containing the relative specific genes used for different mixture types has been attached in bed folders.


mkdir VBSE/VB_1000G/sample-ex
mkdir VBSE/SE_1000G/sample-ex
bed_file="VBSE/VB_specific.bed" 
vcf_file="sample-ex.EAS.dpmaf.recode.vcf" 
#extract the SNP from target genes 
if [[ ! -f "$bed_file" ]]; then
    echo "BED file not found: $bed_file"
    exit 1
fi

if [[ ! -f "$vcf_file" ]]; then
    echo "VCF file not found: $vcf_file"
    exit 1
fi

while read -r chr start end region_name; do
vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out VBSE/VB_1000G/sample-ex/"${region_name//\//_}"
done < "$bed_file"

bed_file="VBSE/SE_specific.bed"

if [[ ! -f "$bed_file" ]]; then
    echo "BED file not found: $bed_file"
    exit 1
fi

if [[ ! -f "$vcf_file" ]]; then
    echo "VCF file not found: $vcf_file"
    exit 1
fi

while read -r chr start end region_name; do
vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out VBSE/SE_1000G/sample-ex/"${region_name//\//_}"
done < "$bed_file"

#clean the files
find VBSE/VB_1000G/sample-ex -type f -size -14860c | xargs rm -f
find VBSE/SE_1000G/sample-ex -type f -size -14860c | xargs rm -f
#transform the files
ls VBSE/VB_1000G/sample-ex/*.recode.vcf | while read id; do grep -v '^##' $id > VBSE/VB_1000G/sample-ex/$(basename $id ".recode.vcf").translated.vcf; done
ls VBSE/SE_1000G/sample-ex/*.recode.vcf | while read id; do grep -v '^##' $id > VBSE/SE_1000G/sample-ex/$(basename $id ".recode.vcf").translated.vcf; done
rm -rf VBSE/VB_1000G/sample-ex/*.recode.vcf VBSE/SE_1000G/sample-ex/*.recode.vcf
#record the SNP using gene names
ls VBSE/VB_1000G/sample-ex/*.translated.vcf | while read id;
do
filename_prefix=$(basename $id ".translated.vcf")
awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
/^#/ {print; next}
{
    $7 = $3;      
    $3 = prefix;  
    print
}' $id > VBSE/VB_1000G/sample-ex/translated_$(basename $id ".translated.vcf").vcf
done
ls VBSE/SE_1000G/sample-ex/*.translated.vcf | while read id;
do
filename_prefix=$(basename $id ".translated.vcf")
awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
/^#/ {print; next}
{
    $7 = $3;      
    $3 = prefix;  
    print
}' $id > VBSE/SE_1000G/sample-ex/translated_$(basename $id ".translated.vcf").vcf
done

#Combine the SNPs from all genes and generate the phenotype frequency
cat VBSE/VB_1000G/sample-ex/translated_* > VBSE/VB_1000G/sample-ex/alltranslatedVBSNP.vcf
cat VBSE/SE_1000G/sample-ex/translated_* > VBSE/SE_1000G/sample-ex/alltranslatedSESNP.vcf
grep -v "^#" VBSE/VB_1000G/sample-ex/alltranslatedVBSNP.vcf > VBSE/VB_1000G/sample-ex/test.vcf
grep -v "^#" VBSE/SE_1000G/sample-ex/alltranslatedSESNP.vcf > VBSE/SE_1000G/sample-ex/test.vcf
cat head.vcf VBSE/VB_1000G/sample-ex/test.vcf > VBSE/VB_1000G/sample-ex/alltranslatedVBSNP.vcf
cat head.vcf VBSE/SE_1000G/sample-ex/test.vcf > VBSE/SE_1000G/sample-ex/alltranslatedSESNP.vcf
rm -rf VBSE/VB_1000G/sample-ex/translated_* VBSE/SE_1000G/sample-ex/translated_*
rm -rf VBSE/VB_1000G/sample-ex/test.vcf VBSE/SE_1000G/sample-ex/test.vcf
python pheno_freq.py -v VBSE/VB_1000G/sample-ex/alltranslatedVBSNP.vcf  -p EAS.txt > VBSE/result/sample-ex.VB_EAS_phenofreq.txt
python pheno_freq.py -v VBSE/SE_1000G/sample-ex/alltranslatedSESNP.vcf  -p EAS.txt > VBSE/result/sample-ex.SE_EAS_phenofreq.txt
