#This code is used for the deconvolution of VB-SE mixtures, but could be extended to other multi-body-fluid-mixtures.
#The folders name could be revised and the folders(BF_mixture) and some of codes could be added according to the deconvoluted mixture types.
#The bed files containing the relative specific genes used for different mixture types has been attached in bed folders.


mkdir VBSE/VB_mixture/sample-ex
mkdir VBSE/SE_mixture/sample-ex
bed_file="VBSE/VB_specific.bed" 
vcf_file="sample-ex.dpmaf.recode.vcf"
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
vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out VBSE/VB_mixture/sample-ex/"${region_name//\//_}"
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
vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out VBSE/SE_mixture/sample-ex/"${region_name//\//_}"
done < "$bed_file"

find VBSE/VB_mixture/sample-ex -type f -size -173927c | xargs rm -f
find VBSE/SE_mixture/sample-ex -type f -size -173927c | xargs rm -f
#transform the files
ls VBSE/VB_mixture/sample-ex/*.recode.vcf | while read id; do bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' $id > VBSE/VB_mixture/sample-ex/$(basename $id ".recode.vcf").vcf; done
ls VBSE/SE_mixture/sample-ex/*.recode.vcf | while read id; do bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' $id > VBSE/SE_mixture/sample-ex/$(basename $id ".recode.vcf").vcf; done
rm -f VBSE/VB_mixture/sample-ex/*.recode.vcf VBSE/SE_mixture/sample-ex/*.recode.vcf
#clean the files
ls VBSE/VB_mixture/sample-ex/*.vcf | while read id;
do
cut -f 1-8 $id > VBSE/VB_mixture/sample-ex/middle.vcf
sed 's/\//|/g' VBSE/VB_mixture/sample-ex/middle.vcf > VBSE/VB_mixture/sample-ex/middle2.vcf
cat head_for_mix.vcf VBSE/VB_mixture/sample-ex/middle2.vcf > VBSE/VB_mixture/sample-ex/$(basename $id ".vcf").translated.vcf
done
ls VBSE/SE_mixture/sample-ex/*.vcf | while read id;
do
cut -f 1-8 $id > VBSE/SE_mixture/sample-ex/middle.vcf
sed 's/\//|/g' VBSE/SE_mixture/sample-ex/middle.vcf > VBSE/SE_mixture/sample-ex/middle2.vcf
cat head_for_mix.vcf VBSE/SE_mixture/sample-ex/middle2.vcf > VBSE/SE_mixture/sample-ex/$(basename $id ".vcf").translated.vcf
done
rm -rf VBSE/VB_mixture/sample-ex/middle* VBSE/SE_mixture/sample-ex/middle*

#record the SNP using gene names
ls VBSE/VB_mixture/sample-ex/*.translated.vcf | while read id;
do
filename_prefix=$(basename $id ".translated.vcf")
awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
/^#/ {print; next}
{
    $7 = $3;      
    $3 = prefix;  
    print
}' $id > VBSE/VB_mixture/sample-ex/translated_$(basename $id ".translated.vcf").vcf
done
ls VBSE/SE_mixture/sample-ex/*.translated.vcf | while read id;
do
filename_prefix=$(basename $id ".translated.vcf")
awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
/^#/ {print; next}
{
    $7 = $3;      
    $3 = prefix;  
    print
}' $id > VBSE/SE_mixture/sample-ex/translated_$(basename $id ".translated.vcf").vcf
done

#Combine the SNPs from all genes and generate the phenotype
cat VBSE/VB_mixture/sample-ex/translated_* > VBSE/VB_mixture/sample-ex/allVBmix.vcf
cat VBSE/SE_mixture/sample-ex/translated_* > VBSE/SE_mixture/sample-ex/allSEmix.vcf
grep -v "^#" VBSE/VB_mixture/sample-ex/allVBmix.vcf > VBSE/VB_mixture/sample-ex/middle.vcf
grep -v "^#" VBSE/SE_mixture/sample-ex/allSEmix.vcf > VBSE/SE_mixture/sample-ex/middle.vcf
cat head_for_mix.vcf VBSE/VB_mixture/sample-ex/middle.vcf > VBSE/VB_mixture/sample-ex/allVBmix.vcf
cat head_for_mix.vcf VBSE/SE_mixture/sample-ex/middle.vcf > VBSE/SE_mixture/sample-ex/allSEmix.vcf
rm -rf VBSE/VB_mixture/sample-ex/translated_* VBSE/SE_mixture/sample-ex/translated_*
rm -rf VBSE/VB_mixture/sample-ex/middle* VBSE/SE_mixture/sample-ex/middle*
python pheno_freq_for_mix.py -v VBSE/VB_mixture/sample-ex/allVBmix.vcf -p mixture.txt > VBSE/result/sample-ex.VBmix_pheno.txt
python pheno_freq_for_mix.py -v VBSE/SE_mixture/sample-ex/allSEmix.vcf -p mixture.txt > VBSE/result/sample-ex.SEmix_pheno.txt



