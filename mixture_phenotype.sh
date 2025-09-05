#This code is used for the deconvolution of two-body-fluid mixtures, but could be extended to other multi-body-fluid-mixtures.


#!/bin/sh

# SCRIPT FUNCTION:
# Accepts a sample ID (-s), two additional IDs (-cF, -cS), and two BED file paths (-c, -d),
# and processes these inputs.

# --- PART 1: CORE PROCESSING LOGIC ---
# This function now receives all five parameters for processing.
mix_phenotype() {
  # Assign the five incoming parameters to local variables for readability and ease of use
  local id="$1"
  local bf_Fst="$2"
  local bf_Sec="$3"
  local bed_file_Fst="$4"
  local bed_file_Sec="$5"

    
 #
 # v v v v v v v v v v v v v v v v v v v v v v v v v v v
    mkdir mixture/result
    mkdir mixture/${bf_Fst}/${sample_id}
    mkdir mixture/${bf_Sec}/${sample_id}
    bed_file="${bed_file_Fst}" 
    vcf_file="${sample_id}.dpmaf.recode.vcf"
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
    vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out mixture/${bf_Fst}/${sample_id}/"${region_name//\//_}"
    done < "$bed_file"
    
    bed_file="${bed_file_Sec}"
    if [[ ! -f "$bed_file" ]]; then
        echo "BED file not found: $bed_file"
        exit 1
    fi
    
    if [[ ! -f "$vcf_file" ]]; then
        echo "VCF file not found: $vcf_file"
        exit 1
    fi
    
    while read -r chr start end region_name; do
    vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out mixture/${bf_Sec}/${sample_id}/"${region_name//\//_}"
    done < "$bed_file"
    
    find mixture/${bf_Fst}/${sample_id} -type f -size -173927c | xargs rm -f
    find mixture/${bf_Sec}/${sample_id} -type f -size -173927c | xargs rm -f
    #transform the files
    ls mixture/${bf_Fst}/${sample_id}/*.recode.vcf | while read id; do bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' $id > mixture/${bf_Fst}/${sample_id}/$(basename $id ".recode.vcf").vcf; done
    ls mixture/${bf_Sec}/${sample_id}/*.recode.vcf | while read id; do bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT]\n' $id > mixture/${bf_Sec}/${sample_id}/$(basename $id ".recode.vcf").vcf; done
    rm -f mixture/${bf_Fst}/${sample_id}/*.recode.vcf mixture/${bf_Sec}/${sample_id}/*.recode.vcf
    #clean the files
    ls mixture/${bf_Fst}/${sample_id}/*.vcf | while read id;
    do
    cut -f 1-8 $id > mixture/${bf_Fst}/${sample_id}/middle.vcf
    sed 's/\//|/g' mixture/${bf_Fst}/${sample_id}/middle.vcf > mixture/${bf_Fst}/${sample_id}/middle2.vcf
    cat head_for_mix.vcf mixture/${bf_Fst}/${sample_id}/middle2.vcf > mixture/${bf_Fst}/${sample_id}/$(basename $id ".vcf").translated.vcf
    done
    ls mixture/${bf_Sec}/${sample_id}/*.vcf | while read id;
    do
    cut -f 1-8 $id > mixture/${bf_Sec}/${sample_id}/middle.vcf
    sed 's/\//|/g' mixture/${bf_Sec}/${sample_id}/middle.vcf > mixture/${bf_Sec}/${sample_id}/middle2.vcf
    cat head_for_mix.vcf mixture/${bf_Sec}/${sample_id}/middle2.vcf > mixture/${bf_Sec}/${sample_id}/$(basename $id ".vcf").translated.vcf
    done
    rm -rf mixture/${bf_Fst}/${sample_id}/middle* mixture/${bf_Sec}/${sample_id}/middle*
    
    #record the SNP using gene names
    ls mixture/${bf_Fst}/${sample_id}/*.translated.vcf | while read id;
    do
    filename_prefix=$(basename $id ".translated.vcf")
    awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
    /^#/ {print; next}
    {
        $7 = $3;      
        $3 = prefix;  
        print
    }' $id > mixture/${bf_Fst}/${sample_id}/translated_$(basename $id ".translated.vcf").vcf
    done
    ls mixture/${bf_Sec}/${sample_id}/*.translated.vcf | while read id;
    do
    filename_prefix=$(basename $id ".translated.vcf")
    awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
    /^#/ {print; next}
    {
        $7 = $3;      
        $3 = prefix;  
        print
    }' $id > mixture/${bf_Sec}/${sample_id}/translated_$(basename $id ".translated.vcf").vcf
    done
    
    #Combine the SNPs from all genes and generate the phenotype
    cat mixture/${bf_Fst}/${sample_id}/translated_* > mixture/${bf_Fst}/${sample_id}/${bf_Fst}.allSNP.vcf
    cat mixture/${bf_Sec}/${sample_id}/translated_* > mixture/${bf_Sec}/${sample_id}/${bf_Sec}.allSNP.vcf
    grep -v "^#" mixture/${bf_Fst}/${sample_id}/${bf_Fst}.allSNP.vcf > mixture/${bf_Fst}/${sample_id}/middle.vcf
    grep -v "^#" mixture/${bf_Sec}/${sample_id}/${bf_Sec}.allSNP.vcf > mixture/${bf_Sec}/${sample_id}/middle.vcf
    cat head_for_mix.vcf mixture/${bf_Fst}/${sample_id}/middle.vcf > mixture/${bf_Fst}/${sample_id}/${bf_Fst}.allSNP.vcf
    cat head_for_mix.vcf mixture/${bf_Sec}/${sample_id}/middle.vcf > mixture/${bf_Sec}/${sample_id}/${bf_Sec}.allSNP.vcf
    rm -rf mixture/${bf_Fst}/${sample_id}/translated_* mixture/${bf_Sec}/${sample_id}/translated_*
    rm -rf mixture/${bf_Fst}/${sample_id}/middle* mixture/${bf_Sec}/${sample_id}/middle*
    python pheno_freq_for_mix.py -v mixture/${bf_Fst}/${sample_id}/${bf_Fst}.allSNP.vcf -p mixture.txt > mixture/result/${sample_id}.${bf_Fst}.mix_pheno.txt
    python pheno_freq_for_mix.py -v mixture/${bf_Sec}/${sample_id}/${bf_Sec}.allSNP.vcf -p mixture.txt > mixture/result/${sample_id}.${bf_Sec}.mix_pheno.txt
 # ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
 #

  echo "Processing complete."
}


# --- PART 2: COMMAND-LINE ARGUMENT PARSING ---
# Initialize all required variables as empty strings
SAMPLE_ID=""
BF_FST=""
BF_SEC=""
BED_FST=""
BED_SEC=""

# Manually parse command-line arguments to support multi-character options
while [ "$#" -gt 0 ]; do
  # Ensure the option is followed by a value
  if [ -z "$2" ]; then
    echo "Error: Option '$1' requires an argument." >&2
    exit 1
  fi

  case "$1" in
    -s)  SAMPLE_ID="$2"; shift 2 ;;
    -cF) BF_FST="$2"; shift 2 ;;     
    -cS) BF_SEC="$2"; shift 2 ;;    
    -bF)  BED_FST="$2"; shift 2 ;;
    -bS)  BED_SEC="$2"; shift 2 ;;
    *)
      echo "Error: Unknown or invalid option: $1" >&2
      echo "Usage: $0 -s <SampleID> -cF <ID_A> -cS <ID_B> -c <path/to/C.bed> -d <path/to/D.bed>" >&2
      exit 1
      ;;
  esac
done


# --- PART 3: MAIN EXECUTION LOGIC ---
# Check if all required parameters have been provided
if [ -z "$SAMPLE_ID" ] || [ -z "$ID_A" ] || [ -z "$ID_B" ] || [ -z "$BED_C" ] || [ -z "$BED_D" ]; then
  echo "Error: Missing required arguments." >&2
  echo "Usage: $0 -s <SampleID> -cF <ID_A> -cS <ID_B> -c <path/to/C.bed> -d <path/to/D.bed>" >&2
  exit 1
else
  # If all parameters have been provided, call the processing function
  mix_phenotype "$SAMPLE_ID" "$ID_A" "$ID_B" "$BED_C" "$BED_D"
fi

echo "----------------------------------------"
echo "==> Script execution finished."



