Fraginae Transcriptome Pipeline Key Steps

### Table of content
#### 1. Assembly
#### 2. Filtering and reduce redundancy  
#### 3. Assessment
#### 4. Quantification and DEG
#### 5. Annotation
#### 6. OrthoFinder and Gene List
#### 7. GO enrichment
#### 8. Candidate Genes

### 1. Assembly

```bash
##### This script is for meta-transcriptome de novo assembly on CU Research Computing Facilities using Trinity.
##### Ruiqi Li
##### Activate python 3 and Trinity environments before running
##### source /curc/sw/anaconda3/latest
##### conda activate Trinity
##### Run script
##### sbatch scriptname

#SBATCH --nodes=1
#SBATCH --time=7-00:00
#SBATCH --partition=smem
#SBATCH --ntasks=48
#SBATCH --job-name=suezensis_ruiqi
#SBATCH --output=SU_trinity.%j.out

source /curc/sw/anaconda3/latest
conda activate Trinity

Trinity --seqType fq --trimmomatic --max_memory 1024G --bflyHeapSpaceMax 1024G --bflyCalculateCPU --left RawSequencingfiles1_1.fq.gz,RawSequencingfiles2_1.fq.gz,...,RawSequencingfilesN_1.fq.gz --right RawSequencingfiles1_2.fq.gz,RawSequencingfiles2_2.fq.gz,...,RawSequencingfilesN_2.fq.gz --CPU 48


### clean heads for trinity.fasta
cat trinity_out_dir/trinity.fasta | sed 's/\s.*$//' > new_fasta_file_headers_trimmedclean.fasta

```

### 2. Filtering and reduce redundancy

```bash
#Prepare local Fragum anmd Symbiodiniaceae Database
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Fragum_whitleyi/latest_assembly_versions/GCA_948146395.1_xbFraWhit1.1/GCA_948146395.1_xbFraWhit1.1_genomic.fna.gz
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Fragum_fragum/latest_assembly_versions/GCA_946902895.1_xbFraFrag2.1/GCA_946902895.1_xbFraFrag2.1_genomic.fna.gz
sudo bash -c 'cat GCA_946902895.1_xbFraFrag2.1_genomic.fna GCA_948146395.1_xbFraWhit1.1_genomic.fna > Fragum_ref.fna'


# Download symbiodiniaceae genomes
# Breviolum minutum: GCA_000507305.1 ASM50730v1 scaffolds: 21,899 contigs: 33,816 N50: 58,535 L50: 3,117
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/507/305/GCA_000507305.1_ASM50730v1/GCA_000507305.1_ASM50730v1_genomic.fna.gz
# Symbiodinium pilosum: 	GCA_905231905.1 Spil_CCMP2461 scaffolds: 48,302 contigs: 142,969 N50: 17,505 L50: 17,600
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/231/905/GCA_905231905.1_Spil_CCMP2461/GCA_905231905.1_Spil_CCMP2461_genomic.fna.gz
#Symbiodinium necroappetens: 	GCA_905231915.1 Snec_CCMP2469 scaffolds: 104,583 contigs: 157,685 N50: 11,421 L50: 18,401
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/231/915/GCA_905231915.1_Snec_CCMP2469/GCA_905231915.1_Snec_CCMP2469_genomic.fna.gz
# Symbiodinium natans: 	GCA_905221605.1 Snat_CCMP2548 scaffolds: 2,855 contigs: 4,262 N50: 358,021 L50: 639
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/221/605/GCA_905221605.1_Snat_CCMP2548/GCA_905221605.1_Snat_CCMP2548_genomic.fna.gz
# Symbiodinium_microadriaticum: GCA_001939145.1_ASM193914v1
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Symbiodinium_microadriaticum/latest_assembly_versions/GCA_001939145.1_ASM193914v1/GCA_001939145.1_ASM193914v1_genomic.fna.gz
# Symbiodinium kawagutii: GCA_009767595.1_ASM976759v1
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/767/595/GCA_009767595.1_ASM976759v1/GCA_009767595.1_ASM976759v1_genomic.fna.gz
# Cladocopium goreaui: GCA_947184155.1_Cgoreaui_SCF055-01
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/184/155/GCA_947184155.1_Cgoreaui_SCF055-01/GCA_947184155.1_Cgoreaui_SCF055-01_genomic.fna.gz
# Cladocopium sp: 	GCA_003297045.1 SymC ver 1.0 scaffolds: 6,576 contigs: 87,131 N50: 38,666 L50: 5,029
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/297/045/GCA_003297045.1_SymC_ver_1.0/GCA_003297045.1_SymC_ver_1.0_genomic.fna.gz
# Symbiodinium sp. clade A Y106 (dinoflagellates): GCA_003297005.1_SymA_ver_1.0_genomic.fna.gz
sudo wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/297/005/GCA_003297005.1_SymA_ver_1.0/GCA_003297005.1_SymA_ver_1.0_genomic.fna.gz

### combine files
sudo bash -c 'cat GCA_000507305.1_ASM50730v1_genomic.fna GCA_003297045.1_SymC_ver_1.0_genomic.fna GCA_905231905.1_Spil_CCMP2461_genomic.fna GCA_001939145.1_ASM193914v1_genomic.fna GCA_009767595.1_ASM976759v1_genomic.fna GCA_905231915.1_Snec_CCMP2469_genomic.fna GCA_003297005.1_SymA_ver_1.0_genomic.fna  GCA_905221605.1_Snat_CCMP2548_genomic.fna GCA_947184155.1_Cgoreaui_SCF055-01_genomic.fna > ../SYMB.fna'
##### add string to the headers
sudo bash -c "sed 's/>/>SYMB_/' SYMB.fna > SYMB_renamed.fna"
#### combine Fragum and symbionts
sudo bash -c 'cat Fragum_ref.fna SYMB_renamed.fna > Fragum_Symb_ref.fna'

# Align the holobiont transcriptome to the local database
nohup mmseqs easy-search --search-type 3 --threads 34 $Query $db SU_alnRes.m8 tmp&
# get the top hit and exclude symbiont hit (SYMB_)
nohup sort -k 1,1 -u SU_alnRes.m8 | grep -v "SYMB_" | cut -f1 > SU_Fragum_tophit_IDs&

# retrieve fasta sequences from the indexed transcriptomic assembly
cdbfasta $Query
cat SU_Fragum_tophit_IDs | cdbyank ${Query}.cidx > 01_host_SU.fna
```

```bash
# 1) remove rRNA by bbmap;
mapPacBio.sh \
in=$rawfile \
ref=$rnaref \
outm=rRNA_${species}.fna \
outu=mRNA_host_${species}_pre.fna

cdbfasta $rawfile

grep ">" mRNA_host_${species}_pre.fna | sed s'/>//' | cut -d ' ' -f1 | sort -u | cdbyank ${rawfile}.cidx > 02_mRNA_host_${species}.fna

# 2) find Open Reading Frames and keep the longest;
TransDecoder.LongOrfs -t 02_mRNA_host_${species}.fna

blastp -query 02_mRNA_host_${species}.fna.transdecoder_dir/longest_orfs.pep  \
    -db $uniprotdb  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 35 > ${species}_blastp.outfmt6

hmmscan --cpu 35 -E 1e-10 --domtblout ${species}_pfam.domtblout $pfamdb 02_mRNA_host_${species}.fna.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t 02_mRNA_host_${species}.fna --retain_pfam_hits ${species}_pfam.domtblout --retain_blastp_hits ${species}_blastp.outfmt6
mv 02_mRNA_host_${species}.fna.transdecoder.cds 03_mRNA_host_${species}.fna.transdecoder.cds

# 3) CD-Hit 90% threshold;
cd-hit-est \
-i 03_mRNA_host_${species}.fna.transdecoder.cds \
-o 04_mRNA_host_${species}.fna.transdecoder.cdhit90.cds \
-c 0.90 \
-n 10 \
-M 0 \
-T 35

```

### 3. Assessment

```bash
#busco
busco -m transcriptome -i $rawfile -o 01_${species}_host_mollusca  -l mollusca_odb10 -c 33
busco -m transcriptome -i $rawfile -o 01_${species}_host_metazoa  -l metazoa_odb10 -c 33

busco -m transcriptome -i 02_mRNA_host_${species}.fna -o 02_mRNA_${species}_mollusca  -l mollusca_odb10 -c 33
busco -m transcriptome -i 02_mRNA_host_${species}.fna -o 02_mRNA_${species}_metazoa  -l metazoa_odb10 -c 33

busco -m transcriptome -i 03_mRNA_host_${species}.fna.transdecoder.cds -o 03_${species}_transdecoder_mollusca  -l mollusca_odb10 -c 33
busco -m transcriptome -i 03_mRNA_host_${species}.fna.transdecoder.cds -o 03_${species}_transdecoder_metazoa  -l metazoa_odb10 -c 33

busco -m transcriptome -i 04_mRNA_host_${species}.fna.transdecoder.cdhit90.cds -o 04_${species}_transdecoder_cdhit90_mollusca  -l mollusca_odb10 -c 33
busco -m transcriptome -i 04_mRNA_host_${species}.fna.transdecoder.cdhit90.cds -o 04_${species}_transdecoder_cdhit90_metazoa  -l metazoa_odb10 -c 33

busco -m transcriptome -i 05_unigene_${species}.fasta -o 05_unigene_${species}_mollusca  -l mollusca_odb10 -c 33
busco -m transcriptome -i 05_unigene_${species}.fasta -o 05_unigene_${species}_metazoa  -l metazoa_odb10 -c 33

mkdir ${species}_BUSCO_mollusca_summaries
cp *_mollusca/short_summary.specific.*.txt ./${species}_BUSCO_mollusca_summaries/
python3 ~/miniconda3/pkgs/busco-4.1.3-py_0/python-scripts/generate_plot.py -wd ${species}_BUSCO_mollusca_summaries
mv ${species}_BUSCO_mollusca_summaries/busco_figure.png ${species}_BUSCO_mollusca_summaries/${species}_BUSCO_mollusca_figure.png

# abyss-fac
abyss-fac 02_mRNA_host_${species}.fna >> ${species}_abyss.txt
abyss-fac 03_mRNA_host_${species}.fna.transdecoder.cds >> ${species}_abyss.txt
abyss-fac 04_mRNA_host_${species}.fna.transdecoder.cdhit90.cds >> ${species}_abyss.txt
abyss-fac 05_unigene_${species}.fasta >> ${species}_abyss.txt
```

### 4. Quantification and DEG
```bash

# 1.Use align_and_estimate_abundance.pl to quantify unigenes
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/align_and_estimate_abundance.pl \
--transcripts unigene_${species}_final.fasta \
--seqType fq \
--samples_file ${species}_samples.txt \
--est_method kallisto \
--prep_reference \
--kallisto_add_opts "-t 34"
# 2.Building Gene Expression Matrices
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/bin/abundance_estimates_to_matrix.pl --est_method kallisto \
--gene_trans_map none \
--out_prefix ${species}_all \
--name_sample_by_basedir \
--quant_files ${species}_abundance_files.list
### Counting Numbers of Expressed Transcripts or Genes
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
${species}_all.isoform.TPM.not_cross_norm | tee ${species}_all.isoform_matrix.TPM.not_cross_norm.counts_by_min_TPM
# Use any TPM>5 as threshold to filter the samples
# Retained 26558 / 75300 = 35.27% of total transcripts.
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/util/filter_low_expr_transcripts.pl --matrix FS_all.isoform.TPM.not_cross_norm --transcripts unigene_FS_final.fasta --min_expr_any 5 > TPM5_unigene_${species}.fasta

# 3.Quality Check Samples and Biological Replicates in all sample dataset
## run PtR
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/Analysis/DifferentialExpression/PtR \
--matrix ${species}_all.isoform.counts.matrix \
--samples ${species}_Quant_samples.txt --log2 --CPM \
--min_rowSums 5 \
--compare_replicates \

## 4.2.1 heatmap
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/Analysis/DifferentialExpression/PtR \
--matrix ${species}_all.isoform.counts.matrix \
--min_rowSums 5 \
-s ${species}_Quant_samples.txt --log2 --CPM --sample_cor_matrix

###PCA
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/Analysis/DifferentialExpression/PtR \
--matrix ${species}_all.isoform.counts.matrix \
-s ${species}_Quant_samples.txt \
--min_rowSums 5 --log2 \
--CPM --center_rows \
--prin_comp 3

```

```r
# Plot heatmap, PCA figures
#Create matrix DESEQ2
Transcrip_all_counts = as.data.frame(lapply(read.table("TG_all.isoform.counts.matrix", row.names = 1, header = T), as.integer))
str(Transcrip_all_counts)
# Define replicates for both conditions to be compared
condition <- factor(read.table("TG_QC_samples.txt", sep = "\t", header= FALSE) %>% select(V1) %>% pull(V1))

# define order
#condition <- factor( condition, levels = c("AM", "AF", "BM", "BF", "CM", "CF"))
# Create a countDataSet
dds = DESeqDataSetFromMatrix(Transcrip_all_counts, DataFrame(condition), design = ~condition)
# Pre-filtering
# keep only rows that have at least 5 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# Estimate normalization factors
dds = estimateSizeFactors(dds)
sizeFactors(dds)

# Estimate dispersion
dds = estimateDispersions(dds)

# Differential analysis
dds  = nbinomWaldTest(dds)

rld = vst(dds, blind=TRUE)

pcaData <- plotPCA(rld, intgroup = c("condition"),  returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$Tissue <- substr(pcaData$group, 2, 2)
pcaData$Treatment <- substr(pcaData$group, 1, 1)
pcaData <- pcaData %>%
  mutate(across('Tissue', str_replace, 'F', 'Foot')) %>%
  mutate(across('Tissue', str_replace, 'M', 'Mantle')) %>%
  mutate(across('Treatment', str_replace, 'A', 'Ambient')) %>%
  mutate(across('Treatment', str_replace, 'B', 'Dim')) %>%
  mutate(across('Treatment', str_replace, 'C', 'Dark'))
pcaData$Treatment <- factor(pcaData$Treatment, levels = c("Ambient", "Dim", "Dark"))

pca<-ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Tissue)) +
  geom_point(size=2, position=position_jitter(height=2, width=2), alpha=0.8) + theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(title = "Principle Component Analysis") +
  coord_fixed(ratio = 1.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text = element_text(face = "bold", size=12), axis.title = element_text(face = "bold", size=13),
        legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 12)) +
  scale_color_manual(values = c("#D55E00", "#E69F00", "#999999"))


# Sample Correlation Heatmap
vsd <- vst(dds, blind=TRUE)
sampleDists <- dist(t(assay(vsd)))

samples<- read.table("TG_QC_samples.txt", sep = "\t", header= FALSE) %>% select(V2) %>% pull(V2)
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samples
colnames(sampleDistMatrix) <- samples
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sample_corr<- pheatmap(sampleDistMatrix,
         #clustering_distance_rows=sampleDists,
         #clustering_distance_cols=sampleDists,
         col=colors)
```

```bash
# DEG
## isoform list (to keep):
genelist=expr5_unigene_${species}.fasta
#get the gene list (min_expr_any>3)
grep ">" $genelist | sed 's/\s.*$//' | sed 's/>//' > expr5_unigene_${species}.geneID
#get the headers
head -1 $count >expr5_${species}_all.isoform.counts.matrix
# get the count matrix on the isoform list
grep -wFf expr5_unigene_${species}.geneID $count >> expr5_${species}_all.isoform.counts.matrix
## Step 1: Identifying DEs
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/bin/run_DE_analysis.pl \
--matrix expr5_${species}_all.isoform.counts.matrix \
--samples_file ${species}_DE_samples_noCF.txt \
--contrasts ${species}_contrasts.txt \
--output . \
--method DESeq2
#Extracting and Clustering DEs
/home/ruiqi/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/bin/analyze_diff_expr.pl \
--matrix expr5_${species}_all.isoform.TMM.EXPR.matrix \
--samples ${species}_DE_samples_noCF.txt \
--output ${species}_P${Pvalue}_C${C} \
--max_genes_clust 30000 \
--P 5e-2 \
--C 1
```


```r
#Plot DE numbers
library(tidyverse)
library(ggplot2)
library(grid)

df <- read.table(file="DE_bargraph.csv", sep = ",", header = TRUE)
# Reorder Comparisons
df$Comparison <- factor(df$Comparison, levels = c("AMvsCM", "BMvsCM", "AMvsBM"))
p <- ggplot(data=df, aes(x=Comparison, y=n, fill=Species)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y=element_blank())+
  labs(y = "Number of DEGs")+
  scale_fill_manual(values=c('#90EE90','#B0E2FF'))+
  #coord_flip()+
  geom_text(aes(x=Comparison, y=n, label=abs(n)), position = position_dodge(width = 1), vjust = ifelse(df$n >= 0, -0.5, 1.5))+
  scale_y_continuous(labels=abs)
```
### 5.Annotation
```bash
emapper.py -m diamond --output_dir $outputdir --translate -o ${species} --cpu 15 -i $input
```

### 6. OrthoFinder and Gene List

```bash
# Run OrthoFinder; -d for amino acids
nohup orthofinder -d -f ./&
```

```r
Ortho <- read.table("Orthogroups.tsv", sep = "\t", header=TRUE)
## split to 2 species
FS <- Ortho %>% pivot_longer(
  cols=starts_with("unigene"),
  values_to = "Gene") %>% filter(name == "unigene_FS_final")
## split multiple genes in one orthogroup
FS_tidy<-separate_rows(FS, Gene, sep=", ")
## Reorder columns: Gene  Orthogroup; delete empty values
FS_table <- FS_tidy %>% select(Gene, Orthogroup) %>% mutate_all(na_if,"") %>% na.omit
## check duplicates
FS_dupes <- FS_table %>% get_dupes(Gene)
write.table(FS_table, file = "FS_GeneOrthoTable.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
```

Genelist to OrthoList, then to Venn Diagram
```bash
for file in FS*.genelist
do
  name="${file%.genelist}"
  echo $name
  grep -Fwf $file FS_GeneOrthoTable.txt | cut -f2 | sort -u > ${name}.ortholist
  rm $file
done
```

```bash
#DE genelist to list 1,2,3
# List 1: Up-regulated genes in FS: AM_CM, excluding TG: AM_CM: 269-2=267
grep -vFwf TG_AM_vs_CM_UP.ortholist FS_AM_vs_CM_UP.ortholist | wc -l
grep -vFwf TG_AM_vs_CM_UP.ortholist FS_AM_vs_CM_UP.ortholist > ../Symbiosis_lists/List1_OG_IDs

# List 2: Up-regulated genes in SU: AM_AF, excluding TG: AM_AF: 1665 - 130 =1535
grep -vFwf TG_AM_vs_AF_UP.ortholist FS_AM_vs_AF_UP.ortholist | wc -l
grep -vFwf TG_AM_vs_AF_UP.ortholist FS_AM_vs_AF_UP.ortholist > ../Symbiosis_lists/List2_OG_IDs

# List 3 Intersection of List 1 and 2
grep -Fwf List1_OG_IDs List2_OG_IDs > List3_OG_IDs

# Step 2:
# List 1: Up-regulated genes in SU: AM_CM, excluding TG: AM_CM

awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' ../Orthogroup_lists/FS_GeneOrthoTable.txt ../Venn_UpSet_figures/FS_AM_vs_CM_UP.genelist | grep -wFf List1_OG_IDs | awk '{print $1}' > list1.genelist

# List 2: Up-regulated genes in SU: AM_AF, excluding TG: AM_AF
awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' ../Orthogroup_lists/FS_GeneOrthoTable.txt ../Venn_UpSet_figures/FS_AM_vs_AF_UP.genelist | grep -wFf List2_OG_IDs | awk '{print $1}' > list2.genelist
# List 3: use the list3 OG ID to get genelists from AM_vs_CM-AM-UP genelist and AM_vs_AF-AM-UP genelist, then get the intersection
awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' ../Orthogroup_lists/FS_GeneOrthoTable.txt ../Venn_UpSet_figures/FS_AM_vs_CM_UP.genelist | grep -wFf List3_OG_IDs | awk '{print $1}' > AM_vs_CM_list3.genelist

awk 'NR==FNR{ a[$1]=$2; next }$1 in a{ $2=a[$1]; print }' ../Orthogroup_lists/FS_GeneOrthoTable.txt ../Venn_UpSet_figures/FS_AM_vs_AF_UP.genelist | grep -wFf List3_OG_IDs | awk '{print $1}' > AM_vs_AF_list3.genelist
# get the intersection
grep -Fwf AM_vs_CM_list3.genelist AM_vs_AF_list3.genelist > list3.genelist

#List 4 Novel Genes
rep -Fwf ../Venn_UpSet_figures/FS_AM_vs_AF_UP.genelist FS_Unassigned_IDs > FS_Unassigned_AM_vs_AF_UP
wc -l FS_Unassigned_AM_vs_AF_UP Okay # 1164/3818
# Get novel genes in AM vs CM up
grep -Fwf ../Venn_UpSet_figures/FS_AM_vs_CM_UP.genelist FS_Unassigned_IDs > FS_Unassigned_AM_vs_CM_UP
wc -l FS_Unassigned_AM_vs_CM_UP # 253/590
# intersection
grep -Fwf FS_Unassigned_AM_vs_CM_UP FS_Unassigned_AM_vs_AF_UP > list4.genelist

```

```r
#Plot the boxplots
library(tidyverse)
############################################# Data Tidying #############################################
# Read in genelist
list3 <- read.table(file="list3.genelist_filtered", sep = "\t", col.names = c("FS_Gene"))
# Read in Gene Expression Matrix TMM
FS_expr<-read.table(file = "expr5_FS_all.isoform.TMM.EXPR.matrix", sep="\t",header = TRUE)
colnames(FS_expr)[1] <- "FS_Gene"

### Add Preferred gene names
genename_raw <- read.table(file = "FS_ID2genename", sep="\t",header = FALSE)
colnames(genename_raw)<-c("FS_Gene", "genename")
#  empty lines to NAs, then dsicard lines with NAs
genename <- genename_raw %>% mutate_all(na_if,"") %>% drop_na
# genename-SU_gene_id_TG and SU expression table
list3_genename <-left_join(list3, genename, by='FS_Gene') %>% drop_na
######### In case duplicated gene names #########
list3_genename$genename <- make.unique(list3_genename$genename, sep = ".")

# Combine blast_Sig with expression matrix
list3_expr <- left_join(list3_genename, FS_expr, by='FS_Gene')

# Data tidying
list3_tidy <- list3_expr %>%
  pivot_longer(!c(FS_Gene, genename), names_to = "Sample", values_to = "TMM")
# Treatment=4th character in TG_A_01M or SU_A_aF
list3_tidy$Treatment <- as.factor(substr(list3_tidy$Sample, 4, 4))
# Tissue= last character in TG_A_01M or SU_A_aF
list3_tidy$Tissue <- as.factor(str_sub(list3_tidy$Sample, -1, -1))
# Species= first 2 characters in G_A_01M or SU_A_aF
list3_tidy$Species <- as.factor(str_sub(list3_tidy$Sample, 1, 2))

# Make it formal
list3_tidy <- list3_tidy%>%
  mutate(across('Tissue', str_replace, 'F', 'Foot')) %>%
  mutate(across('Tissue', str_replace, 'M', 'Mantle')) %>%
  mutate(across('Treatment', str_replace, 'A', 'Ambient')) %>%
  mutate(across('Treatment', str_replace, 'B', 'Dim')) %>%
  mutate(across('Treatment', str_replace, 'C', 'Dark'))
list3_tidy$Treatment <- factor(list3_tidy$Treatment, levels=c('Dark', 'Dim', 'Ambient'))
list3_tidy$Tissue <- as.factor(list3_tidy$Tissue)


############################################# Plotting #############################################
###########  single plot ##########
my_colors <- c("#444444", "#E69F00", "#D55E00")

list3_tidy %>% filter(FS_Gene == "FS_gene_25544") %>%
  ggplot(aes(y=TMM, x=Treatment, fill=Treatment)) +
  geom_boxplot() +
  ggtitle("FS_gene_25544")+
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  facet_wrap(.~Tissue, scales="free") +
  #spacing between two plots
  theme(panel.spacing = unit(2, "lines"))
```

### 7. GO enrichment
####7.1 GO_MWU
```bash
cp *.DE_results DE_statistics/
## Then get genename and logFC
for file in *.DE_results
do
  echo $file
  sp=${file:6:2}
  echo $sp
  comp=${file:35:8}
  echo $comp
  # cut the first line, and get genename (1) and logFC (7)
  cat $file | awk 'NR>1' | cut -f1,7 | tr '\t' ',' > ${sp}_${comp}_lfc.csv
  #rm $file
done

For edgeR Results
file=expr5_FS_all.isoform.counts.matrix.AF_vs_CF.edgeR.DE_results
file=expr5_FS_all.isoform.counts.matrix.CM_vs_CF.edgeR.DE_results
echo $file
sp=${file:6:2}
echo $sp
comp=${file:35:8}
echo $comp
# cut the first line, and get genename (1) and logFC (4)
cat $file | awk 'NR>1' | cut -f1,4 | tr '\t' ',' > ${sp}_${comp}_lfc.csv
```

```r
goDatabase="go.obo"
source("gomwu.functions.R")
###########
# SU
###########
goAnnotations="min_expr_3_SU_gene2go.tab"
# Input
go=c("MF", "CC", "BP")
input_files=c("SU_AM_vs_CM_lfc.csv","SU_AM_vs_BM_lfc.csv","SU_BM_vs_CM_lfc.csv","SU_AM_vs_AF_lfc.csv","SU_CM_vs_CF_lfc.csv","SU_AF_vs_CF_lfc.csv")
# Calculating stats
for (goDivision in go){
  for (input in input_files){
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl",
               largest=0.1,
               smallest=5,
               clusterCutHeight=0.25,
    )
  }
}
```
####7.2 TopGO adn Revigo

```r
geneID2GO <- readMappings(file = "FS_gene2go_topGO.tab")
str(head(geneID2GO))

##############################
### FS - Mantle
##############################

# Load the list of interesting genes
SU_AM_vs_CM_UP_geneID <- read.table(file = "FS_AM_vs_AF_UP_genelist", sep="\t")
SU_AM_vs_CM_UP_geneID <- SU_AM_vs_CM_UP_geneID[,1]
# assign 1 for the DEG, 0 for others
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% SU_AM_vs_CM_UP_geneID))   
names(geneList) <- geneNames
str(geneList)
#construct a topGOdata object
GOdata <- new("topGOdata", description="SU_AM_vs_CM_UP BP", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
# GO enrichment test with fisher classic
SU_AM_vs_CM_UP_resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# visualize p-value vs. number of GOs
pvalFis <- score(SU_AM_vs_CM_UP_resultFisher)
head(pvalFis)
hist(pvalFis, 50, xlab = "p-values")
# Basic information on input data
geneData(SU_AM_vs_CM_UP_resultFisher)
# Summarising the results
allRes <- GenTable(GOdata, classic = SU_AM_vs_CM_UP_resultFisher, ranksOf = "classic", topNodes = 1000)
write.table(allRes, "FS_AM_vs_AF_UP_TopGO.txt", sep = "\t", quote = FALSE)

```

```bash
cat expr5_TG_all.isoform.counts.matrix.AM_vs_AF.DESeq2.DE_results.P5e-2_C1.AF-UP.subset | awk 'NR>1' | cut -f1 > TG_AM_vs_AF_DOWN_genelist
cat FS_AM_vs_AF_UP_TopGO.txt | awk 'NR>1' | cut -f2,7 > Revigo/FS_AM_vs_AF_UP_GO2pvalue.txt
```
Website: [Revigo](http://revigo.irb.hr/)

Parameter:
| Parameter     | value     |
| :------------- | :------------- |
|   similarity     |   0.5     |
| Rank by       | pvalue       |
| remove obsolete GO terms       | Yes   |
| What species       | Whole Uniprot Database       |
| semantic similarity measure       | SimRel       |


#### 8. Candidate Genes

```bash
cat Photosymbiosis_genes.txt | cut -f6,7 | tr -d '"' | sed -e '1d'  | tr "\t" "\n" > gene_quote.fasta
conda activate mmseqs2
#pwd: /home/ruiqi/ruiqi_data/Candidate_Gene_Sequence_Transcriptome
# makedb
makeblastdb -in expr5_unigene_FS.fasta -dbtype nucl -parse_seqids -out expr5_unigene_FS
makeblastdb -in expr5_unigene_TG.fasta -dbtype nucl -parse_seqids -out expr5_unigene_TG
## tblastn
tblastn -query gene_quote_Oct2023.fasta \
    -db expr5_unigene_FS \
    -outfmt 6 -evalue 1e-5 -num_threads 30 > FS_blastp_Oct2023.outfmt6

tblastn -query gene_quote_Oct2023.fasta \
    -db expr5_unigene_TG \
    -outfmt 6 -evalue 1e-5 -num_threads 30 > TG_blastp_Oct2023.outfmt6
```
```r
# plot line plots
#  Read in data
FS_all <- read.table("FS_all_Oct2023.tab", sep = "\t", header = TRUE)
FS_all$AM_vs_BM_Sig <- as.factor(FS_all$AM_vs_BM_Sig)
FS_all$AM_vs_CM_Sig <- as.factor(FS_all$AM_vs_CM_Sig)
FS_all$BM_vs_CM_Sig <- as.factor(FS_all$BM_vs_CM_Sig)
FS_all$AF_vs_CF_Sig <- as.factor(FS_all$AF_vs_CF_Sig)
FS_all$Treatment <- factor(FS_all$Treatment, levels=c('Dark', 'Dim', 'Ambient'))
FS_all$Tissue <- as.factor(FS_all$Tissue)
FS_all$Species <- as.factor(FS_all$Species)
# seperate the M F data
FS_Mantle <- FS_all %>% filter(Tissue == "Mantle") %>% select(-c("AF_vs_CF_Sig")) %>%  na.omit()
FS_Foot <- FS_all %>% filter(Tissue == "Foot") %>% select(-c("AM_vs_BM_Sig","AM_vs_CM_Sig","BM_vs_CM_Sig")) %>% na.omit()

# TG
TG_all <- read.table("TG_all_Oct2023.tab", sep = "\t", header = TRUE)
TG_all$AM_vs_BM_Sig <- as.factor(TG_all$AM_vs_BM_Sig)
TG_all$AM_vs_CM_Sig <- as.factor(TG_all$AM_vs_CM_Sig)
TG_all$BM_vs_CM_Sig <- as.factor(TG_all$BM_vs_CM_Sig)
TG_all$AF_vs_CF_Sig <- as.factor(TG_all$AF_vs_CF_Sig)
TG_all$AF_vs_BF_Sig <- as.factor(TG_all$AF_vs_BF_Sig)
TG_all$BF_vs_CF_Sig <- as.factor(TG_all$BF_vs_CF_Sig)
TG_all$Treatment <- factor(TG_all$Treatment, levels=c('Dark', 'Dim', 'Ambient'))
TG_all$Tissue <- as.factor(TG_all$Tissue)
TG_all$Species <- as.factor(TG_all$Species)

# seperate the M F data
TG_Mantle <- TG_all %>% filter(Tissue == "Mantle") %>% select(-c("AF_vs_CF_Sig", "AF_vs_BF_Sig", "BF_vs_CF_Sig")) %>%  na.omit()
TG_Foot <- TG_all %>% filter(Tissue == "Foot") %>% select(-c("AM_vs_BM_Sig","AM_vs_CM_Sig","BM_vs_CM_Sig", "AF_vs_BF_Sig", "BF_vs_CF_Sig")) %>% na.omit()

# Mantle - two species
colnames(FS_Mantle)[2] <- "Gene"
colnames(TG_Mantle)[2] <- "Gene"
Mantle_FS_TG <- rbind(FS_Mantle, TG_Mantle)

# Foot - two species
colnames(FS_Foot)[2] <- "Gene"
colnames(TG_Foot)[2] <- "Gene"
Foot_FS_TG <- rbind(FS_Foot, TG_Foot)

# Only Keep A-M Sig
TG_Mantle <- TG_all %>% filter(Tissue == "Mantle") %>% select(-c("AF_vs_CF_Sig", "AF_vs_BF_Sig", "BF_vs_CF_Sig", "AM_vs_BM_Sig", "BM_vs_CM_Sig")) %>%  na.omit()
FS_Mantle <- FS_all %>% filter(Tissue == "Mantle") %>% select(-c("AF_vs_CF_Sig", "AM_vs_BM_Sig", "BM_vs_CM_Sig")) %>%  na.omit()
colnames(FS_Mantle)[2] <- "Gene"
colnames(TG_Mantle)[2] <- "Gene"
Mantle_FS_TG <- rbind(FS_Mantle, TG_Mantle)

colnames(Mantle_FS_TG)[3] <- "A_vs_C_Sig"
colnames(Foot_FS_TG)[3] <- "A_vs_C_Sig"
all_FS_TG <- rbind(Mantle_FS_TG, Foot_FS_TG)


png("blast_Mantle_FS_TG_overview_Oct2023.png", width = 5000, height = 6000)
Mantle_FS_TG %>%
  ggplot(aes(y=logFC, x=Treatment, group=Gene)) +
  geom_line(aes(color=AM_vs_CM_Sig), size=5)+
  scale_color_manual(values=c('Grey','Red')) +
  theme_bw() +
  facet_wrap(.~query+Species, ncol = 8) +
  theme(strip.text = element_text(size = 60))
dev.off()
```
