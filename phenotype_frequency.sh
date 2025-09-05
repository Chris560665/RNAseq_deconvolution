#This code is used for the deconvolution of two-body-fluid mixtures, but could be extended to other multi-body-fluid-mixtures.


#!/bin/sh

# SCRIPT FUNCTION:
# Accepts a sample ID (-s), two additional IDs (-cF, -cS), and two BED file paths (-c, -d),
# and processes these inputs.

# --- PART 1: CORE PROCESSING LOGIC ---
# This function now receives all five parameters for processing.
mix_phenotype() {
  # Assign the five incoming parameters to local variables for readability and ease of use
  local sample_id="$1"
  local bf_Fst="$2"
  local bf_Sec="$3"
  local bed_file_Fst="$4"
  local bed_file_Sec="$5"

    
 #
 # v v v v v v v v v v v v v v v v v v v v v v v v v v v
    mkdir -p 1000G/result
    mkdir -p 1000G/${bf_Fst}/${sample_id}
    mkdir -p 1000G/${bf_Sec}/${sample_id}
    bed_file="${bed_file_Fst}" 
    vcf_file="${sample_id}.EAS.dpmaf.recode.vcf" 

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
    vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out 1000G/${bf_Fst}/${sample_id}/"${region_name//\//_}"
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
    vcftools --vcf "$vcf_file" --chr "$chr" --from-bp "$start" --to-bp "$end" --recode --recode-INFO-all --out 1000G/${bf_Sec}/${sample_id}/"${region_name//\//_}"
    done < "$bed_file"
    

    #transform the files
    ls 1000G/${bf_Fst}/${sample_id}/*.recode.vcf | while read id; do grep -v '^##' $id > 1000G/${bf_Fst}/${sample_id}/$(basename $id ".recode.vcf").translated.vcf; done
    ls 1000G/${bf_Sec}/${sample_id}/*.recode.vcf | while read id; do grep -v '^##' $id > 1000G/${bf_Sec}/${sample_id}/$(basename $id ".recode.vcf").translated.vcf; done
    rm -rf 1000G/${bf_Fst}/${sample_id}/*.recode.vcf 1000G/${bf_Sec}/${sample_id}/*.recode.vcf
    #record the SNP using gene names
    ls 1000G/${bf_Fst}/${sample_id}/*.translated.vcf | while read id;
    do
    filename_prefix=$(basename $id ".translated.vcf")
    awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
    /^#/ {print; next}
    {
        $7 = $3;      
        $3 = prefix;  
        print
    }' $id > 1000G/${bf_Fst}/${sample_id}/translated_$(basename $id ".translated.vcf").vcf
    done
    ls 1000G/${bf_Sec}/${sample_id}/*.translated.vcf | while read id;
    do
    filename_prefix=$(basename $id ".translated.vcf")
    awk -v prefix="$filename_prefix" 'BEGIN {OFS="\t"}
    /^#/ {print; next}
    {
        $7 = $3;      
        $3 = prefix;  
        print
    }' $id > 1000G/${bf_Sec}/${sample_id}/translated_$(basename $id ".translated.vcf").vcf
    done
    
    #Combine the SNPs from all genes and generate the phenotype frequency
    cat 1000G/${bf_Fst}/${sample_id}/translated_* > 1000G/${bf_Fst}/${sample_id}/alltranslated${bf_Fst}SNP.vcf
    cat 1000G/${bf_Sec}/${sample_id}/translated_* > 1000G/${bf_Sec}/${sample_id}/alltranslated${bf_Sec}SNP.vcf
    grep -v "^#" 1000G/${bf_Fst}/${sample_id}/alltranslated${bf_Fst}SNP.vcf > 1000G/${bf_Fst}/${sample_id}/test.vcf
    grep -v "^#" 1000G/${bf_Sec}/${sample_id}/alltranslated${bf_Sec}SNP.vcf > 1000G/${bf_Sec}/${sample_id}/test.vcf
    cat head.vcf 1000G/${bf_Fst}/${sample_id}/test.vcf > 1000G/${bf_Fst}/${sample_id}/alltranslated${bf_Fst}SNP.vcf
    cat head.vcf 1000G/${bf_Sec}/${sample_id}/test.vcf > 1000G/${bf_Sec}/${sample_id}/alltranslated${bf_Sec}SNP.vcf
    rm -rf 1000G/${bf_Fst}/${sample_id}/translated_* 1000G/${bf_Sec}/${sample_id}/translated_*
    rm -rf 1000G/${bf_Fst}/${sample_id}/test.vcf 1000G/${bf_Sec}/${sample_id}/test.vcf
    python pheno_freq.py -v 1000G/${bf_Fst}/${sample_id}/alltranslated${bf_Fst}SNP.vcf  -p EAS.txt > 1000G/result/${sample_id}.${bf_Fst}_EAS_phenofreq.txt
    python pheno_freq.py -v 1000G/${bf_Sec}/${sample_id}/alltranslated${bf_Sec}SNP.vcf  -p EAS.txt > 1000G/result/${sample_id}.${bf_Sec}_EAS_phenofreq.txt
    
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
      echo "Usage: $0 -s <SampleID> -cF <BF_FST> -cS <BF_SEC> -bF <path/to/BED_FST.bed> -bS <path/to/BED_SEC.bed>" >&2
      exit 1
      ;;
  esac
done


# --- PART 3: MAIN EXECUTION LOGIC ---
# Check if all required parameters have been provided
if [ -z "$SAMPLE_ID" ] || [ -z "$BF_FST" ] || [ -z "$BF_SEC" ] || [ -z "$BED_FST" ] || [ -z "$BED_SEC" ]; then
  echo "Error: Missing required arguments." >&2
  echo "Usage: $0 -s <SampleID> -cF <BF_FST> -cS <BF_SEC> -bF <path/to/BED_FST.bed> -bS <path/to/BED_SEC.bed>" >&2
  exit 1
else
  # If all parameters have been provided, call the processing function
  mix_phenotype "$SAMPLE_ID" "$BF_FST" "$BF_SEC" "$BED_FST" "$BED_SEC"
fi

echo "----------------------------------------"
echo "==> Script execution finished."