# Code used to process the target captured libraries from trimmed fastq to a summary coverage file. 

# Indexing Reference
bwa index -a bwtsw ../AngiomammalsBacREF_ENT0.9.fasta

# Mapping, sorting, retaining mapped reads only, and cleaning directory
while read name;
do

bwa aln -n 3 ../AngiomammalsBacREF_ENT0.9.fasta ../TC_TRIMMED/$name.trimmed.merged.fastq > $name.TC.3max.AlnENTNM.sai
bwa samse ../AngiomammalsBacREF_ENT0.9.fasta $name.TC.3max.AlnENTNM.sai ../TC_TRIMMED/$name.trimmed.merged.fastq > $name.TC.3max.AlnENTNM.bam
samtools sort $name.TC.3max.AlnENTNM.bam -o $name.TC.3max.AlnENTNM.sorted.bam
samtools index $name.TC.3max.AlnENTNM.sorted.bam
samtools flagstat $name.TC.3max.AlnENTNM.sorted.bam > $name.TC.3max.AlnENTNM.sorted.bam.flagstat.txt
samtools view -F 4 -b $name.TC.3max.AlnENTNM.sorted.bam > $name.TC.3max.AlnENTNM.sorted.mapped.bam
samtools index $name.TC.3max.AlnENTNM.sorted.mapped.bam

rm $name.TC.3max.AlnENTNM.sai
rm $name.TC.3max.AlnENTNM.bam
rm $name.TC.3max.AlnENTNM.sorted.bam

done < ../GSSnames.txt

# Getting coverate statistics
while read name;
do
ls *.sorted.mapped.bam > mergedbamlist
/redser4/software/src/BBTools/bbmap/pileup.sh in=$name.TC.3max.AlnENTNM.sorted.mapped.bam out=$name.TC.3max.AlnENTNM.FEB15_2023.cov.txt ref=../AngiomammalsBacREF_ENT0.9.fasta binsize=10000 bincov=$name.merged_.bins.txt overwrite=true
done < ../GSSnames.txt

# We must remove underscores from the coverage files
while read name;
do
tr ' ' '_' < $name.TC.3max.AlnENTNM.cov.txt > $name.TC.3max.AlnENTNM.noundrscrs.cov.txt
done < ../GSSnames.txt

# Then we must remove zeroes
while read name;
do
awk '{ if ($6!=0) print $0}' $name.TC.3max.AlnENTNM.noundrscrs.cov.txt > $name.TC.3max.AlnENTNM.nozeroes.cov.txt
done < ../GSSnames.txt

#And now we compile coverage statistics into a single file (CompileCoverage.R available in this repository)
Rscript CompileCoverage.R

# Versions of software used in manuscript
# BWA | v0.7.17 | Heng Li, Richard Durbin, Fast and accurate short read alignment with Burrows–Wheeler transform, Bioinformatics, Volume 25, Issue 14, July 2009, Pages 1754–1760, https://doi.org/10.1093/bioinformatics/btp324
# Samtools | v1.13 | Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., et al. (2009) The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics, 25:2078-2079. [PMID: 19505943]
# BBTools | v38.00 | BBMap - Bushnell B. - sourceforge.net/projects/bbmap/
# R | v4.2.1 | R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
