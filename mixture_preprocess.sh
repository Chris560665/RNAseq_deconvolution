#This code is used for preprocessing of gzvcf files from mixture samples
#use -s for sample id input, for which sample should be a gzvcf file and input sample id should be the basename of gzvcf file
#The files "All.chr.EAS.recode.vcf.gz" is accessed from the 1000G project and obtained by merging vcf files from all chromosomes and extracting peoples from East Asian populations using vcftools
#The file RNAedit_pos.txt contained the positions of SNPs with editing events accessed from REDIportal.
#The file "imprintgene_pos.txt" contained the position of imprint gene.
#!/bin/sh


# Core Processing Logic ---
mixpreprocess() {
  local id="$1" # Assign the first argument to a local variable 'id'

  # Check if the incoming id is empty
  if [ -z "$id" ]; then
    echo "Error: Processing function did not receive an ID." >&2
    return 1 # Return an error
  fi


  #
  bcftools filter -e  'FORMAT/DP < 60 || FORMAT/GQ < 40' ${id}.vcf.gz -o ${id}.dp60.vcf.gz
  vcftools --gzvcf ${id}.dp60.vcf.gz --exclude-positions RNAedit_pos.txt --recode --out ${id}_RED
  vcftools --vcf ${id}_RED.recode.vcf --exclude-bed imprintgene_pos.txt --recode --out ${id}_REDIMP
  bcftools view ${id}_REDIMP.recode.vcf -O b -o ${id}_REDIMP.recode.vcf.gz
  bcftools index ${id}_REDIMP.recode.vcf.gz -t
  zcat ${id}_REDIMP.recode.vcf.gz | grep -v '^#' | awk '{print $1 "\t" $2}' > ${id}_pos_dp.txt
  sed -i 's/chr//g' ${id}_pos_dp.txt 
  vcftools --gzvcf All.chr.EAS.recode.vcf.gz --positions ${id}_pos_dp.txt --recode --out ${id}.EAS.1000G
  vcftools --vcf ${id}.EAS.1000G.recode.vcf --maf 0.1 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out ${id}.EAS.dpmaf
  grep -v '^#' ${id}.EAS.dpmaf.recode.vcf |awk '{print $1 "\t" $2}' > ${id}_pos_dpmaf.txt
  awk '{$1="chr"$1; print}' ${id}_pos_dpmaf.txt > ${id}_chrpos_dpmaf.txt
  vcftools --gzvcf ${id}_REDIMP.recode.vcf.gz --positions ${id}_chrpos_dpmaf.txt --recode --out ${id}.middle
  sed -i 's/,<NON_REF>//g' ${id}.middle.recode.vcf
  vcftools --vcf ${id}.middle.recode.vcf --max-alleles 2 --recode --out ${id}.dpmaf
  grep -v '^#' ${id}.dpmaf.recode.vcf |awk '{print $1 "\t" $2}' > ${id}_pos_dpmaf.txt
  sed -i 's/chr//g' ${id}_pos_dpmaf.txt
  vcftools --vcf ${id}.EAS.1000G.recode.vcf --positions ${id}_pos_dpmaf.txt --recode --out ${id}.EAS.dpmaf
  #

  echo "processing complete."
}

# Initialize a variable to store the value specified by the -s parameter, id should be a basename of gzvcf file
SPECIFIC_ID=""

# Use getopts to parse the -s option
while getopts "s:" opt; do
  case "$opt" in
    s)
      # Store the argument value of the -s option in the SPECIFIC_ID variable
      SPECIFIC_ID="$OPTARG"
      ;;
    \?)
      # Handle invalid options
      echo "Invalid option: -$OPTARG" >&2
      echo "Usage: $0 -s <specific_id>" >&2
      exit 1
      ;;
  esac
done


# --- Part 3: Main Execution Logic ---
# In this simplified version, the -s parameter is mandatory.

if [ -z "$SPECIFIC_ID" ]; then
  # If the SPECIFIC_ID variable is empty, it means the user did not provide the -s parameter
  echo "Error: The -s parameter must be used to specify an ID." >&2
  echo "Usage: $0 -s <specific_id>" >&2
  exit 1
else
  # If the user provided the -s parameter, call the processing function
  mixpreprocess "$SPECIFIC_ID"
fi

echo "----------------------------------------"
echo "==> Script execution finished."


