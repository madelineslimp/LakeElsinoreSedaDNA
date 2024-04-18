# Trimming script used for shotgun and target capture libraries 

while read name;
do
/redser4/software/fastp -l 30 -q 25 -y -Y 30 --detect_adapter_for_pe -i $name*R1_001.fastq.gz -I $name*R2_001.fastq.gz -o $name.trimmed.R1.fastq.gz -O $name.trimmed.R2.fastq -m --merged_out $name.trimmed.merged.fastq.gz
done < GSSnames.txt

# Versions of software
# FastP | v0.23.1 | Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
