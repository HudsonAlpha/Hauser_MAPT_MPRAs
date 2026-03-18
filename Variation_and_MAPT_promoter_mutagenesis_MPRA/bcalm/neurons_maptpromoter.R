library(tidyverse)
library(ggplot2)
library(devtools)
#devtools::install_github("kircherlab/BCalm")
library(BCalm)


#loadcounts
counts <- read_tsv("/cluster/home/rhauser/MPRAsnakeflow/results/experiments/rmhex03Count/reporter_experiment.barcode.Neurons.fromFile.default.min_oligo_threshold_10.tsv")

#remove counts for regions with no label
counts <- counts %>% filter(oligo_name != "no_BC")

counts <- counts %>% dplyr::rename(name = oligo_name)
counts <- counts %>% dplyr::rename(Barcode = barcode)
#counts <- counts %>% dplyr::rename(dna_count_1 = 'DNA 1 (condition HEKs, replicate 1)')
#counts <- counts %>% dplyr::rename(dna_count_2 = 'DNA 2 (condition HEKs, replicate 2)')
#counts <- counts %>% dplyr::rename(dna_count_3 = 'DNA 3 (condition HEKs, replicate 3)')
#counts <- counts %>% dplyr::rename(dna_count_4 = 'DNA 4 (condition HEKs, replicate 4)')
#counts <- counts %>% dplyr::rename(rna_count_1 = 'RNA 1 (condition HEKs, replicate 1)')
#counts <- counts %>% dplyr::rename(rna_count_2 = 'RNA 2 (condition HEKs, replicate 2)')
#counts <- counts %>% dplyr::rename(rna_count_3 = 'RNA 3 (condition HEKs, replicate 3)')
#counts <- counts %>% dplyr::rename(rna_count_4 = 'RNA 4 (condition HEKs, replicate 4)')

#counts <- counts %>% dplyr::select(barcode, name, dna_count_1, rna_count_1, dna_count_2, rna_count_2, dna_count_3, rna_count_3, dna_count_4, rna_count_4)

#get variant mapping
#filter counts for just mapt
maptcounts <- counts %>% filter(str_detect(name, regex("mapt", ignore_case = TRUE)))
maptcounts <- maptcounts %>% filter(!str_detect(name, regex("scram", ignore_case = TRUE)))

#load key
#key <- read_tsv("/cluster/home/rhauser/files/rmhex03/design/variant_reference_key.tsv")
#key$ALT <- key$id
#key <- key %>% dplyr::rename(ID = id)
#key <- key %>% dplyr::rename(REF = control)
#key <- key %>% select(ID, REF, ALT)

#load key
key <- read_tsv("/cluster/home/rhauser/files/rmhex03/design/fiveallelesmpradesign.tsv")

#filter key for just mapt
maptkey <- key %>% filter(str_detect(id, regex("mapt", ignore_case = TRUE)))
maptkey <- maptkey %>% filter(!str_detect(id, regex("scram", ignore_case = TRUE)))
newkey <- maptkey %>%
  group_by(element) %>%
  summarize(
    ID   = dplyr::first(element),
    REF  = id[alleleA],
Last login: Wed Mar 18 12:29:21 on ttys000
(base) rebeccahauser@HAIB1122 ~ % ssh cluster
##################################################################
This is a node only for submitting jobs to the cluster; in order
to ensure high availability, jobs that exceed certain CPU and
memory thresholds are automatically killed.

If you need to transfer files to the cluster, compile software, or
run small-scale tests, please use one of the compile nodes at
compile.cluster.haib.org.
##################################################################

########################################################################
.______        ______     ______  __  ___ ____    ____      ___   
|   _  \      /  __  \   /      ||  |/  / \   \  /   /     / _ \  
|  |_)  |    |  |  |  | |  ,----'|  '  /   \   \/   /     | (_) | 
|      /     |  |  |  | |  |     |    <     \_    _/       \__, | 
|  |\  \----.|  `--'  | |  `----.|  .  \      |  |           / /  
| _| `._____| \______/   \______||__|\__\     |__|          /_/   

########################################################################


Last login: Mon Mar 16 15:51:07 2026 from 172.24.232.152
File system space used: 42TB
    projects/ADFTD 630GB
    home/rhauser 41TB
(base) [rhauser@login02 ~]$ ls
DER-08a_hg38_eQTL.significant.chr17.txt  DNAsumtestdesign2.txt    RNA_Atlas_Register.csv  demux.out       files      go                 meme-5.5.8         modules                                       seq_r_env.yml      testpileup3.pdf
DER-08a_hg38_eQTL.significant.txt        MPRAflow                 Rplots.pdf              encodev4.bb     frombelle  heatmap_200kb.pdf  meme-5.5.8.tar.gz  mprasnakeflow.out                             slurm-8230767.out  wu_0.bam
DNAsumtest.tsv                           MPRAsnakeflow            allsnvs.tsv             fastqlist.txt   genomes    maptcombined.tsv   miniconda3         rmhex02_neurons_3reps_2tail_neglog10pval.bed  testpileup.pdf     wu_0.sorted.bam
DNAsumtestdesign.txt                     MotifbreakerResults.tsv  basespace_dowload.out   fastqlist2.txt  github     meme               miniforge3         scripts                                       testpileup2.pdf
(base) [rhauser@login02 ~]$ cd frombelle/
(base) [rhauser@login02 frombelle]$ ls
for_becky
(base) [rhauser@login02 frombelle]$ cd for_becky/
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ nano mapping_paired_end_reads_to_design_file_2023Feb27.sh 
(base) [rhauser@login02 for_becky]$ cd test_association_custom_script/
(base) [rhauser@login02 test_association_custom_script]$ ls
barcodeTable.txt.gz  designfile_unique.fa.amb.gz  designfile_unique.fa.bwt.gz   designfile_unique.fa.fai.gz  designfile_unique.fa.pac.gz  mapped_reads_1M_2023March13.bam.gz    mapped_reads_1M_2023March13.sam  read2_1Mreads.fastq
designfile.fa.gz     designfile_unique.fa.ann.gz  designfile_unique.fa.dict.gz  designfile_unique.fa.gz      designfile_unique.fa.sa.gz   mapped_reads_1M_2023March13.bedpe.gz  read1_1Mreads.fastq
(base) [rhauser@login02 test_association_custom_script]$ head barcodeTable.txt.gz 
Vh?fbarcodeTable.txt????,??4
???&x'?????ɉA???t???T??ҵT?R?vr3H???????\W????????????????????????=???g?{?{?s?5????z??~???????-??\a?k<???O????x??ϟ???"?/QG???????^Ͽ??y?[Ͼ?????g#????I<??ϟ?7[9?Vs????ρ????D?:G1????o4?^?_?s?￟?T?"?~??<[x%?~?D??*????_{??S??c????>??~7?\?s?F?q`???sϿޯP?#V_???V?[Y?T?V???????Cb-9b???????K??7aҰ߿??*?U?|?W??Z?s~^?????x~??????+?r?x??;X?`???
                                                              ??l?z"?(???6??
K??@?X
???5~?b??{w#(H?????>~_??̳x???`??l??R????P??|?OP?|0??US???Tkq???޼?&X?\???Z?
(6/?Z??{?/?Y,??????߯?zM??d?5?D?lS??N?????֦łC?j?#???
z?܄?M?D???K?	(:?!??~??=?i??T?Mi?c	?ln????????|??d???????L?2·?wjvBѰ?X0M??:???<???frn??>??]???Tܞ3a	>?g?4?W[8???F?????D?F?&^?v?SY?UX??$ޗ|?J?o?#H?H?3??x?	?2???Ӽ5?YK?	?n{~?w=g9ͷ1?w?l?U??$?.?TUP??Q????َy???to?Kx??????>?\$ŧx>?+??n??!?`O?a?9?U?ffL????????tˌ??<?????¨???G?Қ???c|?\H???%??;??7aVۼz?+???6?x?C?wW?<̽,t??(oX??0?%??GlJƦ<?w?Z0??+?Tx???????`mZN?9??k??4
??W??t}>!K?E?(?y??M?\?*?>7?,?t?Ѿ????z+???????XvT?<r?^+??B?W9LG?g????w?b#&?-???K??\?g??G_?^!?#??????︮?%??[??b??w?????l?G??%J????Y??l???
f?????p?????YE?,??XX>?R?K??h????͟???V?ʥ&?EH?lb#?????P,{?_??5??f>BqC=-?yӮ?7?y??O???????Q1T?n?qA?j?g??HLd~-
                                        ?u??8??
??9?gy??ߍޣ?gK??*??
???c??[3X??zlb?in??D?Gh???Q????x??˯?z??7??y?????rܖX????????kw-5*?U?V?y??Z?M]??B??g???x???J?,???A???X?k[>??_?h?*??
(base) [rhauser@login02 test_association_custom_script]$ cd ..
(base) [rhauser@login02 for_becky]$ cd testdesignfastq/
(base) [rhauser@login02 testdesignfastq]$ ls
bamfile.bam  bedfile.paired.bed  bedfile.sorted.uniq.rangeFilter.sizeFilter.bed  designfile.fa  samfile.sam  test_r1.fastq  test_r2.fastq
(base) [rhauser@login02 testdesignfastq]$ cd ..
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ cd testing_header_script/
(base) [rhauser@login02 testing_header_script]$ ls
extract_header_info_forBecky_2022November03.sh  head_HT7GJBGXL_s1.all_info.txt  head_HT7GJBGXL_s1_1_BC_X.fastq  head_HT7GJBGXL_s1_2_BC_X.fastq
(base) [rhauser@login02 testing_header_script]$ cd ..
(base) [rhauser@login02 for_becky]$ cd Assoc_Basic/
(base) [rhauser@login02 Assoc_Basic]$ ls
output
(base) [rhauser@login02 Assoc_Basic]$ cd output/
(base) [rhauser@login02 output]$ ls
assoc_basic
(base) [rhauser@login02 output]$ cd assoc_basic/
(base) [rhauser@login02 assoc_basic]$ ls
assoc_basic_barcode_counts.pickle                                  assoc_basic_filtered_coords_to_barcodes.pickle  assoc_basic_original_counts.png  design_rmIllegalChars.fa    label_rmIllegalChars.txt
assoc_basic_barcodes_per_candidate-no_repeats-no_jackpots.feather  assoc_basic_filtered_counts.png                 count_fastq.txt                  filtered_count_summary.txt  original_count_summary.txt
(base) [rhauser@login02 assoc_basic]$ cd ..
(base) [rhauser@login02 output]$ ls
assoc_basic
(base) [rhauser@login02 output]$ cd ..
(base) [rhauser@login02 Assoc_Basic]$ ls
output
(base) [rhauser@login02 Assoc_Basic]$ cd ..
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ cd barcode_location_assignment_2023April17/
(base) [rhauser@login02 barcode_location_assignment_2023April17]$ ls
bedout_qscore30plus.bed  headerInfo_in_qscore30plus.bed
(base) [rhauser@login02 barcode_location_assignment_2023April17]$ cd ..
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ cd ..
(base) [rhauser@login02 frombelle]$ ls
for_becky
(base) [rhauser@login02 frombelle]$ cd ..
(base) [rhauser@login02 ~]$ ls
DER-08a_hg38_eQTL.significant.chr17.txt  DNAsumtestdesign2.txt    RNA_Atlas_Register.csv  demux.out       files      go                 meme-5.5.8         modules                                       seq_r_env.yml      testpileup3.pdf
DER-08a_hg38_eQTL.significant.txt        MPRAflow                 Rplots.pdf              encodev4.bb     frombelle  heatmap_200kb.pdf  meme-5.5.8.tar.gz  mprasnakeflow.out                             slurm-8230767.out  wu_0.bam
DNAsumtest.tsv                           MPRAsnakeflow            allsnvs.tsv             fastqlist.txt   genomes    maptcombined.tsv   miniconda3         rmhex02_neurons_3reps_2tail_neglog10pval.bed  testpileup.pdf     wu_0.sorted.bam
DNAsumtestdesign.txt                     MotifbreakerResults.tsv  basespace_dowload.out   fastqlist2.txt  github     meme               miniforge3         scripts                                       testpileup2.pdf
(base) [rhauser@login02 ~]$ cd files
(base) [rhauser@login02 files]$ ls
andersonrogers2023supplementaldata  encode  fastqlist.txt  from_henry  from_nick  from_others  graphs  hic_rizzardi  practice  rmhex01  rmhex02  rmhex03  rogers2024supplementaldata  slurm-8238174.out  slurm-8238956.out  temp
(base) [rhauser@login02 files]$ cd practice/
(base) [rhauser@login02 practice]$ ls
(base) [rhauser@login02 practice]$ cd ..
(base) [rhauser@login02 files]$ cd from_others/
(base) [rhauser@login02 from_others]$ ls
ADSP_chr17_biallelic_regions0005.vcf.gz  ADSP_chr17_biallelic_regions20.vcf  ADSP_chr17_biallelic_regions25.vcf  MAPT_DegCre_Trilineage_grant_GWAS.csv  justtheimportantinfo_norm_variants25.tsv
(base) [rhauser@login02 from_others]$ cd ..
(base) [rhauser@login02 files]$ ls
andersonrogers2023supplementaldata  encode  fastqlist.txt  from_henry  from_nick  from_others  graphs  hic_rizzardi  practice  rmhex01  rmhex02  rmhex03  rogers2024supplementaldata  slurm-8238174.out  slurm-8238956.out  temp
(base) [rhauser@login02 files]$ cd temp
(base) [rhauser@login02 temp]$ ls
ADSP_becky_chr17_snps.vcf.gz                                   ADSP_becky_chr17_snps_dbSNP-156_wTOPMed_wAnn_wCADD-1.6.vcf.gz.tbi  rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed
ADSP_becky_chr17_snps.vcf.gz.tbi                               ADSP_becky_chr17_snps_snpeff_stats.genes.txt                       rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed.onlyup
ADSP_becky_chr17_snps_dbSNP-156.vcf.gz                         ADSP_becky_chr17_snps_snpeff_stats.html                            rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed.onlyup.sorted
ADSP_becky_chr17_snps_dbSNP-156.vcf.gz.tbi                     ADSP_becky_chr17_snps_vt-noCHR.vcf.gz                              rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed.onlyup.sorted.merge
ADSP_becky_chr17_snps_dbSNP-156_wTOPMed.vcf.gz                 ADSP_becky_chr17_snps_vt-noCHR.vcf.gz.tbi                          rmhex02_neurons_analysis2_50cutoff_neglog10pval.bed.onlyup.sorted.merge.chr17.bed
ADSP_becky_chr17_snps_dbSNP-156_wTOPMed.vcf.gz.tbi             all.onlyup.chr17.bed                                               slurm-3051374.out
ADSP_becky_chr17_snps_dbSNP-156_wTOPMed_wAnn.vcf.gz            all.onlyup.chr17.bed.sorted                                        slurm-3073536.out
ADSP_becky_chr17_snps_dbSNP-156_wTOPMed_wAnn.vcf.gz.tbi        all.onlyup.chr17.bed.sorted.merge.bed
ADSP_becky_chr17_snps_dbSNP-156_wTOPMed_wAnn_wCADD-1.6.vcf.gz  annotate_vcf_fromjared.sh
(base) [rhauser@login02 temp]$ cd ..
(base) [rhauser@login02 files]$ cd rmhex01
(base) [rhauser@login02 rmhex01]$ ls
250bindesignfile.bed  allrepsfastqs           barcodeassocfastqs                  getfasta_fordesignfile.sh  mpraflowBarcodeAssoc2       rep1fastqs_try2_20230223                             rmhex01_neurons_analysis2_50cutoff_sighits.bed
250bindesignfile.fa   bacslabel.bed           barcodeassocfastqs_files_HT7GJBGXL  makedesignfile             mpraflowBarcodeAssoc3       rep3fastqs                                           sangerseq
MPRAflowCount_2reps   bacslabel_original.bed  dictionary250.pickle                makedesignfile2            nonmergedbarcodeassocalign  rmhex01_heks_analysis2_50cutoff_neglog10pval.bed     slurm-4032189.out
MPRAflowCount_3reps   barcodeallign           forGEO                              makedesignfile3            rep1_reseq_20240201         rmhex01_heks_analysis2_50cutoff_sighits.bed          slurm-4032212.out
MPRAsnakeflowCount    barcodeallignfastmap    gapscheck                           mpraflowBarcodeAssoc       rep1fastqs_20230220         rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed
(base) [rhauser@login02 rmhex01]$ cd nonmergedbarcodeassocalign/
(base) [rhauser@login02 nonmergedbarcodeassocalign]$ ls
newsamfile.sam  samfile.sam
(base) [rhauser@login02 nonmergedbarcodeassocalign]$ cd ..
(base) [rhauser@login02 rmhex01]$ cd barcodeallign
(base) [rhauser@login02 barcodeallign]$ ls
QC          bedcut.bed       bedout.bed     bedsorted_w_QS.bed  bsub.out       fastaout.fa                      filtunmapped_w_qs.bed  forvis        onlyqs40up.bed  range.res.length.count  test.callable.bed  uniqbed_w_QS.bed
bamout.bam  bedcut_w_QS.bed  bedsorted.bed  bedunique.bed       bsubbedqs.txt  fastaoutsizeandrangefiltered.fa  first1000sam.sam       headfasta.fa  onlyqs60.bed    samfile.sam             test.depth.bed
(base) [rhauser@login02 barcodeallign]$ cd ..
(base) [rhauser@login02 rmhex01]$ cd barcodeallignfastmap/
(base) [rhauser@login02 barcodeallignfastmap]$ ls
busublog.txt  bwa_fastmap.matches
(base) [rhauser@login02 barcodeallignfastmap]$ cd ..
(base) [rhauser@login02 rmhex01]$ cd barcodeassocfastqs
(base) [rhauser@login02 barcodeassocfastqs]$ ls
7109-NC-0001.txt                                   HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  dsl_wget.sh                                     read2.fastq.gz    slurm-749848.out
HT7GJBGXL_s1_1_BC_X.fastq                          HT7GJBGXL_s1_2_BC_X.fastq                                    QC                                                           extract_header_info_forBecky_2022November03.sh  slurm-749150.out  slurm-750690.out
HT7GJBGXL_s1_1_BC_X.fastq.gz.md5sum                HT7GJBGXL_s1_2_BC_X.fastq.gz.md5sum                          barcodeinfoextract                                           gzip.bsub.out                                   slurm-749154.out  slurm-750693.out
HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq            barcodeinfoextract2                                          read1.fastq.gz                                  slurm-749845.out
(base) [rhauser@login02 barcodeassocfastqs]$ cd barcodeinfoextract
(base) [rhauser@login02 barcodeinfoextract]$ ls
HT7GJBGXL_s1.all_info.txt       alluniquebarcodes.txt  barcodesread1.fastq  bsubrscript.out     first25ktable.txt  slurm-1062953.out  slurm-1063005.out  slurm-1063222.out  toydata.txt             uniquebccount.py   uniquebccount3.py
HT7GJBGXL_s1.all_info.txt.save  barcodeprepcmd.R       barcodesread2.fastq  first1000table.txt  slurm-1062952.out  slurm-1062963.out  slurm-1063068.out  test               uniquebarcodecount.csv  uniquebccount2.py
(base) [rhauser@login02 barcodeinfoextract]$ cd test/
(base) [rhauser@login02 test]$ ls
barcodesread1.fastq  barcodesread2.fastq  toydata.txt
(base) [rhauser@login02 test]$ cd ..
(base) [rhauser@login02 barcodeinfoextract]$ ls
HT7GJBGXL_s1.all_info.txt       alluniquebarcodes.txt  barcodesread1.fastq  bsubrscript.out     first25ktable.txt  slurm-1062953.out  slurm-1063005.out  slurm-1063222.out  toydata.txt             uniquebccount.py   uniquebccount3.py
HT7GJBGXL_s1.all_info.txt.save  barcodeprepcmd.R       barcodesread2.fastq  first1000table.txt  slurm-1062952.out  slurm-1062963.out  slurm-1063068.out  test               uniquebarcodecount.csv  uniquebccount2.py
(base) [rhauser@login02 barcodeinfoextract]$ cd ..
(base) [rhauser@login02 barcodeassocfastqs]$ ls
7109-NC-0001.txt                                   HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  dsl_wget.sh                                     read2.fastq.gz    slurm-749848.out
HT7GJBGXL_s1_1_BC_X.fastq                          HT7GJBGXL_s1_2_BC_X.fastq                                    QC                                                           extract_header_info_forBecky_2022November03.sh  slurm-749150.out  slurm-750690.out
HT7GJBGXL_s1_1_BC_X.fastq.gz.md5sum                HT7GJBGXL_s1_2_BC_X.fastq.gz.md5sum                          barcodeinfoextract                                           gzip.bsub.out                                   slurm-749154.out  slurm-750693.out
HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq            barcodeinfoextract2                                          read1.fastq.gz                                  slurm-749845.out
(base) [rhauser@login02 barcodeassocfastqs]$ cd barcodeinfoextract2
(base) [rhauser@login02 barcodeinfoextract2]$ ls
barcodeprepcmd2.R  bc.fastq.gz  bsub.bczip.out  bsub.out
(base) [rhauser@login02 barcodeinfoextract2]$ cd ..
(base) [rhauser@login02 barcodeassocfastqs]$ ls
7109-NC-0001.txt                                   HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  dsl_wget.sh                                     read2.fastq.gz    slurm-749848.out
HT7GJBGXL_s1_1_BC_X.fastq                          HT7GJBGXL_s1_2_BC_X.fastq                                    QC                                                           extract_header_info_forBecky_2022November03.sh  slurm-749150.out  slurm-750690.out
HT7GJBGXL_s1_1_BC_X.fastq.gz.md5sum                HT7GJBGXL_s1_2_BC_X.fastq.gz.md5sum                          barcodeinfoextract                                           gzip.bsub.out                                   slurm-749154.out  slurm-750693.out
HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq            barcodeinfoextract2                                          read1.fastq.gz                                  slurm-749845.out
(base) [rhauser@login02 barcodeassocfastqs]$ cd ..
(base) [rhauser@login02 rmhex01]$ ls
250bindesignfile.bed  allrepsfastqs           barcodeassocfastqs                  getfasta_fordesignfile.sh  mpraflowBarcodeAssoc2       rep1fastqs_try2_20230223                             rmhex01_neurons_analysis2_50cutoff_sighits.bed
250bindesignfile.fa   bacslabel.bed           barcodeassocfastqs_files_HT7GJBGXL  makedesignfile             mpraflowBarcodeAssoc3       rep3fastqs                                           sangerseq
MPRAflowCount_2reps   bacslabel_original.bed  dictionary250.pickle                makedesignfile2            nonmergedbarcodeassocalign  rmhex01_heks_analysis2_50cutoff_neglog10pval.bed     slurm-4032189.out
MPRAflowCount_3reps   barcodeallign           forGEO                              makedesignfile3            rep1_reseq_20240201         rmhex01_heks_analysis2_50cutoff_sighits.bed          slurm-4032212.out
MPRAsnakeflowCount    barcodeallignfastmap    gapscheck                           mpraflowBarcodeAssoc       rep1fastqs_20230220         rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed
(base) [rhauser@login02 rmhex01]$ cd barcodeassocfastqs_files_HT7GJBGXL/
(base) [rhauser@login02 barcodeassocfastqs_files_HT7GJBGXL]$ ls
 HT7GJBGXL_s1_1_BC_X.fastq.gz                           HT7GJBGXL_s1_2_BC_X.fastq.gz                           dsl_wget.sh  'files_HT7GJBGXL_20220928_114024 (1).txt'      multiqc_report_haib22NC7109_19712_batch1.html
 HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq.gz   HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq.gz   files.txt     multiqc_data_haib22NC7109_19712_batch1.json
(base) [rhauser@login02 barcodeassocfastqs_files_HT7GJBGXL]$ cd ..
(base) [rhauser@login02 rmhex01]$ ls
250bindesignfile.bed  allrepsfastqs           barcodeassocfastqs                  getfasta_fordesignfile.sh  mpraflowBarcodeAssoc2       rep1fastqs_try2_20230223                             rmhex01_neurons_analysis2_50cutoff_sighits.bed
250bindesignfile.fa   bacslabel.bed           barcodeassocfastqs_files_HT7GJBGXL  makedesignfile             mpraflowBarcodeAssoc3       rep3fastqs                                           sangerseq
MPRAflowCount_2reps   bacslabel_original.bed  dictionary250.pickle                makedesignfile2            nonmergedbarcodeassocalign  rmhex01_heks_analysis2_50cutoff_neglog10pval.bed     slurm-4032189.out
MPRAflowCount_3reps   barcodeallign           forGEO                              makedesignfile3            rep1_reseq_20240201         rmhex01_heks_analysis2_50cutoff_sighits.bed          slurm-4032212.out
MPRAsnakeflowCount    barcodeallignfastmap    gapscheck                           mpraflowBarcodeAssoc       rep1fastqs_20230220         rmhex01_neurons_analysis2_50cutoff_neglog10pval.bed
(base) [rhauser@login02 rmhex01]$ cd barcodeassocfastqs
(base) [rhauser@login02 barcodeassocfastqs]$ ls
7109-NC-0001.txt                                   HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  dsl_wget.sh                                     read2.fastq.gz    slurm-749848.out
HT7GJBGXL_s1_1_BC_X.fastq                          HT7GJBGXL_s1_2_BC_X.fastq                                    QC                                                           extract_header_info_forBecky_2022November03.sh  slurm-749150.out  slurm-750690.out
HT7GJBGXL_s1_1_BC_X.fastq.gz.md5sum                HT7GJBGXL_s1_2_BC_X.fastq.gz.md5sum                          barcodeinfoextract                                           gzip.bsub.out                                   slurm-749154.out  slurm-750693.out
HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq            barcodeinfoextract2                                          read1.fastq.gz                                  slurm-749845.out
(base) [rhauser@login02 barcodeassocfastqs]$ cd barcodeinfoextract
(base) [rhauser@login02 barcodeinfoextract]$ ls
HT7GJBGXL_s1.all_info.txt       alluniquebarcodes.txt  barcodesread1.fastq  bsubrscript.out     first25ktable.txt  slurm-1062953.out  slurm-1063005.out  slurm-1063222.out  toydata.txt             uniquebccount.py   uniquebccount3.py
HT7GJBGXL_s1.all_info.txt.save  barcodeprepcmd.R       barcodesread2.fastq  first1000table.txt  slurm-1062952.out  slurm-1062963.out  slurm-1063068.out  test               uniquebarcodecount.csv  uniquebccount2.py
(base) [rhauser@login02 barcodeinfoextract]$ cd test/
(base) [rhauser@login02 test]$ ls
barcodesread1.fastq  barcodesread2.fastq  toydata.txt
(base) [rhauser@login02 test]$ cd ..
(base) [rhauser@login02 barcodeinfoextract]$ ls
HT7GJBGXL_s1.all_info.txt       alluniquebarcodes.txt  barcodesread1.fastq  bsubrscript.out     first25ktable.txt  slurm-1062953.out  slurm-1063005.out  slurm-1063222.out  toydata.txt             uniquebccount.py   uniquebccount3.py
HT7GJBGXL_s1.all_info.txt.save  barcodeprepcmd.R       barcodesread2.fastq  first1000table.txt  slurm-1062952.out  slurm-1062963.out  slurm-1063068.out  test               uniquebarcodecount.csv  uniquebccount2.py
(base) [rhauser@login02 barcodeinfoextract]$ cd ..
(base) [rhauser@login02 barcodeassocfastqs]$ cd barcodeinfoextract2
(base) [rhauser@login02 barcodeinfoextract2]$ ls
barcodeprepcmd2.R  bc.fastq.gz  bsub.bczip.out  bsub.out
(base) [rhauser@login02 barcodeinfoextract2]$ cd ..
(base) [rhauser@login02 barcodeassocfastqs]$ ls
7109-NC-0001.txt                                   HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq.gz.md5sum  dsl_wget.sh                                     read2.fastq.gz    slurm-749848.out
HT7GJBGXL_s1_1_BC_X.fastq                          HT7GJBGXL_s1_2_BC_X.fastq                                    QC                                                           extract_header_info_forBecky_2022November03.sh  slurm-749150.out  slurm-750690.out
HT7GJBGXL_s1_1_BC_X.fastq.gz.md5sum                HT7GJBGXL_s1_2_BC_X.fastq.gz.md5sum                          barcodeinfoextract                                           gzip.bsub.out                                   slurm-749154.out  slurm-750693.out
HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001.fastq  HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001.fastq            barcodeinfoextract2                                          read1.fastq.gz                                  slurm-749845.out
(base) [rhauser@login02 barcodeassocfastqs]$ cd QC
(base) [rhauser@login02 QC]$ ls
HT7GJBGXL_s1_1_BC_X_fastqc.html  HT7GJBGXL_s1_1_BC_X_screen.txt                           HT7GJBGXL_s1_2_BC_X_fastqc.html  HT7GJBGXL_s1_2_BC_X_screen.txt                          fastqscreenbsublog.txt
HT7GJBGXL_s1_1_BC_X_fastqc.zip   HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001_fastqc.html  HT7GJBGXL_s1_2_BC_X_fastqc.zip   HT7GJBGXL_s1_2_NOID_15-NOID_10_7109-NC-0001_fastqc.zip
HT7GJBGXL_s1_1_BC_X_screen.html  HT7GJBGXL_s1_1_NOID_15-NOID_10_7109-NC-0001_fastqc.zip   HT7GJBGXL_s1_2_BC_X_screen.html  bsublog.txt
(base) [rhauser@login02 QC]$ cd ~
(base) [rhauser@login02 ~]$ ls
DER-08a_hg38_eQTL.significant.chr17.txt  DNAsumtestdesign2.txt    RNA_Atlas_Register.csv  demux.out       files      go                 meme-5.5.8         modules                                       seq_r_env.yml      testpileup3.pdf
DER-08a_hg38_eQTL.significant.txt        MPRAflow                 Rplots.pdf              encodev4.bb     frombelle  heatmap_200kb.pdf  meme-5.5.8.tar.gz  mprasnakeflow.out                             slurm-8230767.out  wu_0.bam
DNAsumtest.tsv                           MPRAsnakeflow            allsnvs.tsv             fastqlist.txt   genomes    maptcombined.tsv   miniconda3         rmhex02_neurons_3reps_2tail_neglog10pval.bed  testpileup.pdf     wu_0.sorted.bam
DNAsumtestdesign.txt                     MotifbreakerResults.tsv  basespace_dowload.out   fastqlist2.txt  github     meme               miniforge3         scripts                                       testpileup2.pdf
(base) [rhauser@login02 ~]$ cd frombelle/
(base) [rhauser@login02 frombelle]$ ls
for_becky
(base) [rhauser@login02 frombelle]$ cd for_becky/
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ cd barcode_location_assignment_2023April17/
(base) [rhauser@login02 barcode_location_assignment_2023April17]$ ls
bedout_qscore30plus.bed  headerInfo_in_qscore30plus.bed
(base) [rhauser@login02 barcode_location_assignment_2023April17]$ cd ..
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ cd Assoc_Basic/
(base) [rhauser@login02 Assoc_Basic]$ ls
output
(base) [rhauser@login02 Assoc_Basic]$ cd output/
(base) [rhauser@login02 output]$ ls
assoc_basic
(base) [rhauser@login02 output]$ cd assoc_basic/
(base) [rhauser@login02 assoc_basic]$ ls
assoc_basic_barcode_counts.pickle                                  assoc_basic_filtered_coords_to_barcodes.pickle  assoc_basic_original_counts.png  design_rmIllegalChars.fa    label_rmIllegalChars.txt
assoc_basic_barcodes_per_candidate-no_repeats-no_jackpots.feather  assoc_basic_filtered_counts.png                 count_fastq.txt                  filtered_count_summary.txt  original_count_summary.txt
(base) [rhauser@login02 assoc_basic]$ cd ..
(base) [rhauser@login02 output]$ ls
assoc_basic
(base) [rhauser@login02 output]$ cd ..
(base) [rhauser@login02 Assoc_Basic]$ ls
output
(base) [rhauser@login02 Assoc_Basic]$ cd ..
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ ls
Assoc_Basic                  HT7GJBGXL_s1_BC_X.fastqun2               designfile.fa      designfile.fa.bwt  designfile.fa.sa                                      mapping_results_2023Feb27_failed.sam  testing_header_script
HT7GJBGXL_s1_BC_X.fastqjoin  barcode_location_assignment_2023April17  designfile.fa.amb  designfile.fa.fai  mapping_paired_end_reads_to_design_file_2023Feb27.sh  test_association_custom_script
HT7GJBGXL_s1_BC_X.fastqun1   bwa_mapping_2023Feb27_logfile.txt        designfile.fa.ann  designfile.fa.pac  mapping_results_2023Feb27.sam                         testdesignfastq
(base) [rhauser@login02 for_becky]$ cd ~
(base) [rhauser@login02 ~]$ cd find . -name "barcode_association_multiple_support.txt"
-bash: cd: too many arguments
(base) [rhauser@login02 ~]$ find . -name "barcode_association_multiple_support.txt"
(base) [rhauser@login02 ~]$ find . -name "barcode_association_multiple_support.txt"
(base) [rhauser@login02 ~]$ find . -name "barcode_association_multiple_support.txt"
./files/rmhex01/associated_files/barcode_association_multiple_support.txt
(base) [rhauser@login02 ~]$ ls
DER-08a_hg38_eQTL.significant.chr17.txt  DNAsumtestdesign2.txt    RNA_Atlas_Register.csv  demux.out       files      go                 meme-5.5.8         modules                                       seq_r_env.yml      testpileup3.pdf
DER-08a_hg38_eQTL.significant.txt        MPRAflow                 Rplots.pdf              encodev4.bb     frombelle  heatmap_200kb.pdf  meme-5.5.8.tar.gz  mprasnakeflow.out                             slurm-8230767.out  wu_0.bam
DNAsumtest.tsv                           MPRAsnakeflow            allsnvs.tsv             fastqlist.txt   genomes    maptcombined.tsv   miniconda3         rmhex02_neurons_3reps_2tail_neglog10pval.bed  testpileup.pdf     wu_0.sorted.bam
DNAsumtestdesign.txt                     MotifbreakerResults.tsv  basespace_dowload.out   fastqlist2.txt  github     meme               miniforge3         scripts                                       testpileup2.pdf
(base) [rhauser@login02 ~]$ cd files
(base) [rhauser@login02 files]$ ls
andersonrogers2023supplementaldata  encode  fastqlist.txt  from_henry  from_nick  from_others  graphs  hic_rizzardi  practice  rmhex01  rmhex02  rmhex03  rogers2024supplementaldata  slurm-8238174.out  slurm-8238956.out  temp
(base) [rhauser@login02 files]$ ls -lh
total 88K
drwxr-xr-x 1 rhauser cochran-lab   0 Apr 26  2024 andersonrogers2023supplementaldata
drwxr-xr-x 1 rhauser cochran-lab   0 Jun 13  2025 encode
-rw-r--r-- 1 rhauser cochran-lab 69K Dec  5 12:28 fastqlist.txt
drwxr-xr-x 1 rhauser cochran-lab   0 Dec 11 09:35 from_henry
drwxr-xr-x 1 rhauser cochran-lab   0 Nov 19 13:05 from_nick
drwxr-xr-x 1 rhauser cochran-lab   0 Aug 12  2025 from_others
drwxr-xr-x 1 rhauser cochran-lab   0 Aug 11  2025 graphs
drwxr-xr-x 1 rhauser cochran-lab   0 Mar  8  2024 hic_rizzardi
drwxr-xr-x 1 rhauser cochran-lab   0 Dec  4 15:45 practice
drwxr-xr-x 1 rhauser cochran-lab   0 Mar 18 13:53 rmhex01
drwxr-xr-x 1 rhauser cochran-lab   0 Mar 17 09:46 rmhex02
drwxr-xr-x 1 rhauser cochran-lab   0 Mar 12 12:50 rmhex03
drwxr-xr-x 1 rhauser cochran-lab   0 Jun 23  2025 rogers2024supplementaldata
-rw-r--r-- 1 rhauser cochran-lab   0 Dec  5 11:54 slurm-8238174.out
-rw-r--r-- 1 rhauser cochran-lab   0 Dec  5 12:27 slurm-8238956.out
drwxr-xr-x 1 rhauser cochran-lab   0 Apr 24  2024 temp
(base) [rhauser@login02 files]$ cd from_others/
(base) [rhauser@login02 from_others]$ ls
ADSP_chr17_biallelic_regions0005.vcf.gz  ADSP_chr17_biallelic_regions20.vcf  ADSP_chr17_biallelic_regions25.vcf  MAPT_DegCre_Trilineage_grant_GWAS.csv  justtheimportantinfo_norm_variants25.tsv
(base) [rhauser@login02 from_others]$ cd ..
(base) [rhauser@login02 files]$ ls
andersonrogers2023supplementaldata  encode  fastqlist.txt  from_henry  from_nick  from_others  graphs  hic_rizzardi  practice  rmhex01  rmhex02  rmhex03  rogers2024supplementaldata  slurm-8238174.out  slurm-8238956.out  temp
(base) [rhauser@login02 files]$ from_jared
-bash: from_jared: command not found
(base) [rhauser@login02 files]$ mkdir from_jared
(base) [rhauser@login02 files]$ rm -r from_jared
(base) [rhauser@login02 files]$ cd ..
(base) [rhauser@login02 ~]$ ls
DER-08a_hg38_eQTL.significant.chr17.txt  DNAsumtestdesign2.txt    RNA_Atlas_Register.csv  demux.out       files      go                 meme-5.5.8         modules                                       seq_r_env.yml      testpileup3.pdf
DER-08a_hg38_eQTL.significant.txt        MPRAflow                 Rplots.pdf              encodev4.bb     frombelle  heatmap_200kb.pdf  meme-5.5.8.tar.gz  mprasnakeflow.out                             slurm-8230767.out  wu_0.bam
DNAsumtest.tsv                           MPRAsnakeflow            allsnvs.tsv             fastqlist.txt   genomes    maptcombined.tsv   miniconda3         rmhex02_neurons_3reps_2tail_neglog10pval.bed  testpileup.pdf     wu_0.sorted.bam
DNAsumtestdesign.txt                     MotifbreakerResults.tsv  basespace_dowload.out   fastqlist2.txt  github     meme               miniforge3         scripts                                       testpileup2.pdf
(base) [rhauser@login02 ~]$ cd files/rmhex02
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd ..
(base) [rhauser@login02 files]$ sftp subftp@sftp-private.ncbi.nlm.nih.gov

This warning banner provides privacy and security notices consistent with 
applicable federal laws, directives, and other federal guidance for accessing 
this Government system, which includes all devices/storage media attached to 
this system. This system is provided for Government-authorized use only. 
Unauthorized or improper use of this system is prohibited and may result in 
disciplinary action and/or civil and criminal penalties. At any time, and for 
any lawful Government purpose, the government may monitor, record, and audit 
your system usage and/or intercept, search and seize any communication or data 
transiting or stored on this system. Therefore, you have no reasonable 
expectation of privacy. Any communication or data transiting or stored on this 
system may be disclosed or used for any lawful Government purpose.
subftp@sftp-private.ncbi.nlm.nih.gov's password: 
Connected to sftp-private.ncbi.nlm.nih.gov.
sftp> cd uploads/rhauser_hudsonalpha.org_rpYoJM6r
sftp> ls
BACMPRA  
sftp> cd BACMPRA/
sftp> ls
BACMPRA_rep1_DNA_Heks_R1_merged.fq.gz                             BACMPRA_rep1_DNA_Heks_R1_sequencing1.fq.gz                        BACMPRA_rep1_DNA_Heks_R1_sequencing2.fq.gz                        BACMPRA_rep1_DNA_Heks_R2_merged_with_faux_fastqs.fq.gz            
BACMPRA_rep1_DNA_Heks_R2_sequencing1.fq.gz                        BACMPRA_rep1_DNA_Heks_R2_sequencing2.fq.gz                        BACMPRA_rep1_DNA_Heks_R3_merged.fq.gz                             BACMPRA_rep1_DNA_Heks_R3_sequencing2.fq.gz                        
BACMPRA_rep1_DNA_Neurons_R1_merged.fq.gz                          BACMPRA_rep1_DNA_Neurons_R1_sequencing1.fq.gz                     BACMPRA_rep1_DNA_Neurons_R1_sequencing2.fq.gz                     BACMPRA_rep1_DNA_Neurons_R2_merged_with_faux_fastqs.fq.gz         
BACMPRA_rep1_DNA_Neurons_R2_sequencing1.fq.gz                     BACMPRA_rep1_DNA_Neurons_R2_sequencing2.fq.gz                     BACMPRA_rep1_DNA_Neurons_R3_merged.fq.gz                          BACMPRA_rep1_DNA_Neurons_R3_sequencing2.fq.gz                     
BACMPRA_rep1_RNA_Heks_R1_merged.fq.gz                             BACMPRA_rep1_RNA_Heks_R1_sequencing1.fq.gz                        BACMPRA_rep1_RNA_Heks_R1_sequencing2.fq.gz                        BACMPRA_rep1_RNA_Heks_R2_merged_with_faux_fastqs.fq.gz            
BACMPRA_rep1_RNA_Heks_R2_sequencing1.fq.gz                        BACMPRA_rep1_RNA_Heks_R2_sequencing2.fq.gz                        BACMPRA_rep1_RNA_Heks_R3_merged.fq.gz                             BACMPRA_rep1_RNA_Heks_R3_sequencing2.fq.gz                        
BACMPRA_rep1_RNA_Neurons_R1_merged.fq.gz                          BACMPRA_rep1_RNA_Neurons_R1_sequencing1.fq.gz                     BACMPRA_rep1_RNA_Neurons_R1_sequencing2.fq.gz                     BACMPRA_rep1_RNA_Neurons_R2_merged_with_faux_fastqs.fq.gz         
BACMPRA_rep1_RNA_Neurons_R2_sequencing1.fq.gz                     BACMPRA_rep1_RNA_Neurons_R2_sequencing2.fq.gz                     BACMPRA_rep1_RNA_Neurons_R3_merged.fq.gz                          BACMPRA_rep1_RNA_Neurons_R3_sequencing2.fq.gz                     
BACMPRA_rep2_DNA_Heks_R1.fq.gz                                    BACMPRA_rep2_DNA_Heks_R2.fq.gz                                    BACMPRA_rep2_DNA_Heks_R3.fq.gz                                    BACMPRA_rep2_DNA_Neurons_R1.fq.gz                                 
BACMPRA_rep2_DNA_Neurons_R2.fq.gz                                 BACMPRA_rep2_DNA_Neurons_R3.fq.gz                                 BACMPRA_rep2_RNA_Heks_R1.fq.gz                                    BACMPRA_rep2_RNA_Heks_R2.fq.gz                                    
BACMPRA_rep2_RNA_Heks_R3.fq.gz                                    BACMPRA_rep2_RNA_Neurons_R1.fq.gz                                 BACMPRA_rep2_RNA_Neurons_R2.fq.gz                                 BACMPRA_rep2_RNA_Neurons_R3.fq.gz                                 
BACMPRA_rep3_DNA_Heks_R1.fq.gz                                    BACMPRA_rep3_DNA_Heks_R2.fq.gz                                    BACMPRA_rep3_DNA_Heks_R3.fq.gz                                    BACMPRA_rep3_DNA_Neurons_R1.fq.gz                                 
BACMPRA_rep3_DNA_Neurons_R2.fq.gz                                 BACMPRA_rep3_DNA_Neurons_R3.fq.gz                                 BACMPRA_rep3_RNA_Heks_R1.fq.gz                                    BACMPRA_rep3_RNA_Heks_R2.fq.gz                                    
BACMPRA_rep3_RNA_Heks_R3.fq.gz                                    BACMPRA_rep3_RNA_Neurons_R1.fq.gz                                 BACMPRA_rep3_RNA_Neurons_R2.fq.gz                                 BACMPRA_rep3_RNA_Neurons_R3.fq.gz                                 
read1.fastq.gz                                                    read2.fastq.gz                                                    
sftp> mv read1.fastq.gz bcassoc_R1.fq.gz
Invalid command.
sftp> rm read
read1.fastq.gz  read2.fastq.gz  
sftp> rm read1.fastq.gz 
Removing /uploads/rhauser_hudsonalpha.org_rpYoJM6r/BACMPRA/read1.fastq.gz
sftp> rm read2.fastq.gz 
Removing /uploads/rhauser_hudsonalpha.org_rpYoJM6r/BACMPRA/read2.fastq.gz
sftp> quit
(base) [rhauser@login02 files]$ cd ~
(base) [rhauser@login02 ~]$ cd files/rmhex02
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd MPRAflowCount_reps123/
(base) [rhauser@login02 MPRAflowCount_reps123]$ ls
MPRAflow_count_3181183.out  MPRAflow_count_3249377.out               MPRAflow_count_nompranalyze_3249378.out  output_mpranalyze                  rmhex02_count_call_shellscript_nompranalyze.sh  work_mpranalyze
MPRAflow_count_3231441.out  MPRAflow_count_nompranalyze_3181196.out  output                                   rmhex02_count_call_shellscript.sh  work
(base) [rhauser@login02 MPRAflowCount_reps123]$ cd output_mpranalyze/
(base) [rhauser@login02 output_mpranalyze]$ ls
HEK293FT  Neurons
(base) [rhauser@login02 output_mpranalyze]$ cd Neurons/
(base) [rhauser@login02 Neurons]$ ls
1                                                                average_bedgraph3772074.out                              rmhex02_3reps_neurons_neglog10pval.bed
2                                                                controls_pval_hist_analysis2.png                         rmhex02_3reps_neurons_neglog10pval.col5bed.bed
3                                                                controls_pval_hist_prefilt_allcontrols.png               rmhex02_3reps_neurons_neglog10pval.col5bed.partition.bed
Neurons_count.csv                                                correctbedfiles                                          rmhex02_3reps_neurons_neglog10pval.col5bed.sorted.bed
Neurons_final_labeled_counts.txt                                 dna_annot.tsv                                            rmhex02_3reps_neurons_neglog10pval.col5bed.sorted.mean.bed
all_pval_hist_analysis2.png                                      dna_counts.tsv                                           rmhex02_3reps_neurons_neglog10pval.partition.bed
all_pval_hist_analysis_uq_prefilt_nonegcontrols.png              dna_counts_padded.csv                                    rmhex02_3reps_neurons_neglog10pval.sorted.bed
all_pval_hist_analysis_uq_prefilt_nonegcontrols_onetail.png      mpranalyze_3252509.out                                   rmhex02_3reps_neurons_neglog10pval.sorted.mean.bed
all_pval_hist_prefilt_allcontrols.png                            mpranalyze_3252511.out                                   rmhex02_3reps_neurons_sighits_up.bed
all_pval_hist_prefilt_allcontrols_oneside.png                    mpranalyze_3315589.out                                   rmhex02_3reps_neurons_sighits_up_.001.bed
all_pval_hist_prefilt_allcontrols_oneside_25cutoff.png           mpranalyze_3315590.out                                   rmhex02_3reps_neurons_sighits_up_.01.bed
all_pval_hist_prefilt_allcontrols_oneside_mean50.png             mpranalyze_3315835.out                                   rmhex02_neurons_3reps_2tail_neglog10pval.bed
alpha_analysis2.tsv                                              mpranalyze_3751059.out                                   rmhex02_neurons_3reps_2tail_neglog10pval.col5bed.bed
alpha_analysis2_oneside.tsv                                      mpranalyze_3756532.out                                   rmhex02_neurons_3reps_2tail_neglog10pval.partition.bed
alpha_analysis2_oneside_25cutoff.tsv                             mpranalyze_3760421.out                                   rmhex02_neurons_3reps_2tail_neglog10pval.sorted.bed
alpha_analysis2_oneside_mean50.tsv                               mpranalyze_4016078.out                                   rmhex02_neurons_3reps_2tail_neglog10pval.sorted.mean.bed
alpha_analysis_uq_prefilt_nonegcontrols.tsv                      mpranalyze_7646666.out                                   rna_annot.tsv
alpha_analysis_uq_prefilt_nonegcontrols_onetail.tsv              mpranalyze_7646668.out                                   rna_counts.tsv
alpha_boxplot_analysis2.png                                      mpranalyze_callscript_prefilt_nonegcontrols.sh           rna_counts_padded.csv
alpha_boxplot_analysis_prefilt_allcontrols.png                   mpranalyze_callscript_prefilt_nonegcontrols_onesided.sh  testEmpirical_alpha_Neurons_results_labeled_analysis_uq_prefilt_nonegcontrols.tsv
alpha_boxplot_analysis_prefilt_allcontrols_oneside.png           mpranalyze_onesided_prefilt_nonegcontrols.R              testEmpirical_alpha_Neurons_results_labeled_analysis_uq_prefilt_nonegcontrols_onetail.tsv
alpha_boxplot_analysis_prefilt_allcontrols_oneside_25cutoff.png  mpranalyze_twosided_prefilt_nonegcontrols.R              testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside.tsv
alpha_boxplot_analysis_prefilt_allcontrols_oneside_mean50.png    mpranalyzecallscript.sh                                  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_25cutoff.tsv
alpha_boxplot_analysis_uq_prefilt_nonegcontrols.png              mpranalyzecallscript_onesided.sh                         testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
alpha_boxplot_analysis_uq_prefilt_nonegcontrols_onetail.png      mpranalyzecallscript_onesided_25cutoff.sh                testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_mean50.tsv
average_bedgraph3663234.out                                      mpranalyzecallscript_onesided_mean50.sh                  testEmpirical_alpha_Neurons_results_prefilt_allcontrols_oneside.tsv
average_bedgraph3663237.out                                      oligoMPRA_onesided.R                                     testEmpirical_alpha_Neurons_results_prefilt_allcontrols_oneside_25cutoff.tsv
average_bedgraph3663239.out                                      oligoMPRA_onesided_25cutoff.R                            testEmpirical_alpha_Neurons_results_prefilt_allcontrols_oneside_mean50.tsv
average_bedgraph3663241.out                                      oligoMPRA_onesided_mean50.R                              testEmpirical_alpha_rmhex02_Neurons_results_labeled_prefilt_allcontrols.tsv
average_bedgraph3663242.out                                      oligoMPRA_twosided.R                                     testEmpirical_alpha_rmhex02_Neurons_results_prefilt_allcontrols.tsv
average_bedgraph3663243.out                                      padding.sh                                               testEmpirical_alpha_rmhex02_neurons_reps123_twosided_results_analysis2.tsv
average_bedgraph3663244.out                                      padding_3252505.out                                      testEmpirical_alpha_rmhex02_neurons_reps123_twosided_results_labeled_analysis2.tsv
(base) [rhauser@login02 Neurons]$ cd ..
(base) [rhauser@login02 output_mpranalyze]$ ls
HEK293FT  Neurons
(base) [rhauser@login02 output_mpranalyze]$ cd HEK293FT/
(base) [rhauser@login02 HEK293FT]$ ls
1                                                            controls_pval_hist_prefilt_allcontrols.png  mpranalyze_7669399.out                                   rmhex02_3reps_heks_sighits_up_.001.bed
2                                                            controls_pval_hist_prefilt_onetail.png      mpranalyze_7669403.out                                   rmhex02_3reps_heks_sighits_up_.01.bed
3                                                            dna_annot.tsv                               mpranalyze_7669404.out                                   rmhex02_4reps_heks_sighits_up.bed
HEK293FT_count.csv                                           dna_counts.tsv                              mpranalyze_7669406.out                                   rna_annot.tsv
HEK293FT_final_labeled_counts.txt                            dna_counts_padded.csv                       mpranalyze_7674873.out                                   rna_counts.tsv
all_pval_hist_analysis2.png                                  mpranalyze_3252510.out                      mpranalyze_7674874.out                                   rna_counts_padded.csv
all_pval_hist_analysis_uq_prefilt_nonegcontrols_onetail.png  mpranalyze_3315592.out                      mpranalyze_7674875.out                                   testEmpirical_alpha_HEK293FT_results_labeled_analysis_uq_prefilt_nonegcontrols_onetail.tsv
all_pval_hist_analysis_uq_prefilt_nonegcontrols_twotail.png  mpranalyze_3752896.out                      mpranalyze_callscript_prefilt_nonegcontrols.sh           testEmpirical_alpha_HEK293FT_results_labeled_analysis_uq_prefilt_nonegcontrols_twotail.tsv
all_pval_hist_prefilt_allcontrols.png                        mpranalyze_3754295.out                      mpranalyze_callscript_prefilt_nonegcontrols_onesided.sh  testEmpirical_alpha_HEK293FT_results_labeled_prefilt_allcontrols.tsv
all_pval_hist_prefilt_allcontrols_mean50.png                 mpranalyze_3760419.out                      mpranalyze_onesided_prefilt_nonegcontrols.R              testEmpirical_alpha_HEK293FT_results_labeled_prefilt_allcontrols_alllabels.tsv
all_pval_hist_prefilt_onetail.png                            mpranalyze_3760420.out                      mpranalyze_twosided_prefilt_nonegcontrols.R              testEmpirical_alpha_HEK293FT_results_labeled_prefilt_allcontrols_mean50.tsv
alpha_analysis2.tsv                                          mpranalyze_3760525.out                      mpranalyzecallscript.sh                                  testEmpirical_alpha_HEK293FT_results_labeled_prefilt_onetail.tsv
alpha_analysis2_mean50.tsv                                   mpranalyze_3760527.out                      mpranalyzecallscript_25cutoff.sh                         testEmpirical_alpha_HEK293FT_results_prefilt_allcontrols.tsv
alpha_analysis_uq_prefilt_nonegcontrols_onetail.tsv          mpranalyze_3761574.out                      mpranalyzecallscript_mean50.sh                           testEmpirical_alpha_HEK293FT_results_prefilt_allcontrols_mean50.tsv
alpha_analysis_uq_prefilt_nonegcontrols_twotail.tsv          mpranalyze_3885146.out                      oligoMPRA_onesided.R                                     testEmpirical_alpha_HEK293FT_results_prefilt_onetail.tsv
alpha_boxplot_analysis2.png                                  mpranalyze_3886132.out                      oligoMPRA_onesided_25cutoff.R                            testEmpirical_alpha_NA_results_labeled_prefilt_onetail.tsv
alpha_boxplot_analysis_prefilt_allcontrols.png               mpranalyze_3886133.out                      oligoMPRA_onesided_mean50.R                              testEmpirical_alpha_NA_results_prefilt_onetail.tsv
alpha_boxplot_analysis_prefilt_allcontrols_mean50.png        mpranalyze_7669110.out                      oligoMPRA_twosided.R                                     testEmpirical_alpha_rmhex02_neurons_reps123_twosided_results_analysis2.tsv
alpha_boxplot_analysis_uq_prefilt_nonegcontrols_onetail.png  mpranalyze_7669113.out                      padding.sh                                               testEmpirical_alpha_rmhex02_neurons_reps123_twosided_results_labeled_analysis2.tsv
alpha_boxplot_analysis_uq_prefilt_nonegcontrols_twotail.png  mpranalyze_7669123.out                      padding_3252506.out
controls_pval_hist_analysis2.png                             mpranalyze_7669136.out                      rmhex02_3reps_heks_neglog10pval.bed
(base) [rhauser@login02 HEK293FT]$ nano oligoMPRA_onesided.R
(base) [rhauser@login02 HEK293FT]$ cd ..
(base) [rhauser@login02 output_mpranalyze]$ ls
HEK293FT  Neurons
(base) [rhauser@login02 output_mpranalyze]$ cd ..
(base) [rhauser@login02 MPRAflowCount_reps123]$ ls
MPRAflow_count_3181183.out  MPRAflow_count_3249377.out               MPRAflow_count_nompranalyze_3249378.out  output_mpranalyze                  rmhex02_count_call_shellscript_nompranalyze.sh  work_mpranalyze
MPRAflow_count_3231441.out  MPRAflow_count_nompranalyze_3181196.out  output                                   rmhex02_count_call_shellscript.sh  work
(base) [rhauser@login02 MPRAflowCount_reps123]$ nano rmhex02_count_call_shellscript_nompranalyze.sh 
(base) [rhauser@login02 MPRAflowCount_reps123]$ nano rmhex02_count_call_shellscript.sh
(base) [rhauser@login02 MPRAflowCount_reps123]$ ls
MPRAflow_count_3181183.out  MPRAflow_count_3249377.out               MPRAflow_count_nompranalyze_3249378.out  output_mpranalyze                  rmhex02_count_call_shellscript_nompranalyze.sh  work_mpranalyze
MPRAflow_count_3231441.out  MPRAflow_count_nompranalyze_3181196.out  output                                   rmhex02_count_call_shellscript.sh  work
(base) [rhauser@login02 MPRAflowCount_reps123]$ cd ..
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd MPRAflowBCassoc
(base) [rhauser@login02 MPRAflowBCassoc]$ ls
mpraflow.error  mpraflow.out  mpraflow_commentoutslurm.error  mpraflow_commentoutslurm.out  output  rmhex02bcassocmpraflow.sh  rmhex02bcassocmpraflow_slurmcommentout.sh  work
(base) [rhauser@login02 MPRAflowBCassoc]$ cd output/
(base) [rhauser@login02 output]$ ls
assoc_basic_RMHex02
(base) [rhauser@login02 output]$ cd assoc_basic_RMHex02/
(base) [rhauser@login02 assoc_basic_RMHex02]$ ls
count_fastq.txt  design_rmIllegalChars.fa  label_rmIllegalChars.txt
(base) [rhauser@login02 assoc_basic_RMHex02]$ cd ..
(base) [rhauser@login02 output]$ ls
assoc_basic_RMHex02
(base) [rhauser@login02 output]$ cd ..
(base) [rhauser@login02 MPRAflowBCassoc]$ ls
mpraflow.error  mpraflow.out  mpraflow_commentoutslurm.error  mpraflow_commentoutslurm.out  output  rmhex02bcassocmpraflow.sh  rmhex02bcassocmpraflow_slurmcommentout.sh  work
(base) [rhauser@login02 MPRAflowBCassoc]$ cd ..
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing
(base) [rhauser@login02 MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing]$ ls
1_count_bc_nolab_new_try3.sh                                                    RMHex02_bcassoc_allreads_barcodes_per_candidate.feather      design_rmIllegalChars2.fa   slurm-1546767.out  slurm-1549113.out
3_PE_merge2.sh                                                                  RMHex02_bcassoc_allreads_coords_to_barcodes.pickle           exploratoryanalysis         slurm-1546769.out  slurm-1549116.out
4_align_BWA_PE2.sh                                                              RMHex02_bcassoc_allreads_filtered_coords_to_barcodes.pickle  filtered_count_summary.txt  slurm-1546776.out  slurm-1772107.out
5_map_element_barcodesnew.sh                                                    RMHex02_bcassoc_allreads_filtered_coords_to_barcodes.tsv     label_rmIllegalChars.txt    slurm-1546778.out  slurm-1774315.out
6.err                                                                           RMHex02_bcassoc_allreads_filtered_counts.png                 labels.txt                  slurm-1546779.out  slurm-1774561.out
6.out                                                                           RMHex02_bcassoc_allreads_original_counts.png                 mpraflowdict_to_tsv.py      slurm-1548139.out  unique_count_mappedelements.sort.txt
6_filter_barcodesnew.sh                                                         all_fastqjoin_merged.count_bam.txt                           mpraflowdict_to_tsv.sh      slurm-1548727.out  unique_count_mappedelements.txt
RMHex02_bcassoc_allreads_barcode_counts.pickle                                  all_fastqjoin_merged.sorted.bam                              original_count_summary.txt  slurm-1549075.out  uniquemappedelements.txt
RMHex02_bcassoc_allreads_barcodes_per_candidate-no_jackpots.feather             all_fastqjoin_merged.sorted.sam                              samfilefacts.err            slurm-1549076.out
RMHex02_bcassoc_allreads_barcodes_per_candidate-no_repeats-no_jackpots.feather  bamtosam3.sh                                                 samfilefacts.out            slurm-1549077.out
RMHex02_bcassoc_allreads_barcodes_per_candidate-no_repeats.feather              count_fastq.txt                                              samfilefacts_newdup.sh      slurm-1549079.out
(base) [rhauser@login02 MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing]$ cd ..
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd ..
(base) [rhauser@login02 files]$ ls
andersonrogers2023supplementaldata  encode  fastqlist.txt  from_henry  from_nick  from_others  graphs  hic_rizzardi  practice  rmhex01  rmhex02  rmhex03  rogers2024supplementaldata  slurm-8238174.out  slurm-8238956.out  temp
(base) [rhauser@login02 files]$ cd rmhex03
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd bcassoc
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd bcassoc2018
-bash: cd: bcassoc2018: No such file or directory
(base) [rhauser@login02 bcassoc]$ cd bcassoc218
(base) [rhauser@login02 bcassoc218]$ ls
MPRAflow_assoc_218_7293354.out  MPRAflow_assoc_218_called_shorttime2_7566671.out  bcassoc218.sh                            bcassoc218_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_218_7294464.out  MPRAflow_assoc_218_called_shorttime2_7614909.out  bcassoc218_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc218]$ nano bcassoc218.sh
(base) [rhauser@login02 bcassoc218]$ cd ..
(base) [rhauser@login02 bcassoc]$ cd ..
(base) [rhauser@login02 rmhex03]$ cd ..
(base) [rhauser@login02 files]$ cd rmhex02
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing 
(base) [rhauser@login02 MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing]$ ls
1_count_bc_nolab_new_try3.sh                                                    RMHex02_bcassoc_allreads_barcodes_per_candidate.feather      design_rmIllegalChars2.fa   slurm-1546767.out  slurm-1549113.out
3_PE_merge2.sh                                                                  RMHex02_bcassoc_allreads_coords_to_barcodes.pickle           exploratoryanalysis         slurm-1546769.out  slurm-1549116.out
4_align_BWA_PE2.sh                                                              RMHex02_bcassoc_allreads_filtered_coords_to_barcodes.pickle  filtered_count_summary.txt  slurm-1546776.out  slurm-1772107.out
5_map_element_barcodesnew.sh                                                    RMHex02_bcassoc_allreads_filtered_coords_to_barcodes.tsv     label_rmIllegalChars.txt    slurm-1546778.out  slurm-1774315.out
6.err                                                                           RMHex02_bcassoc_allreads_filtered_counts.png                 labels.txt                  slurm-1546779.out  slurm-1774561.out
6.out                                                                           RMHex02_bcassoc_allreads_original_counts.png                 mpraflowdict_to_tsv.py      slurm-1548139.out  unique_count_mappedelements.sort.txt
6_filter_barcodesnew.sh                                                         all_fastqjoin_merged.count_bam.txt                           mpraflowdict_to_tsv.sh      slurm-1548727.out  unique_count_mappedelements.txt
RMHex02_bcassoc_allreads_barcode_counts.pickle                                  all_fastqjoin_merged.sorted.bam                              original_count_summary.txt  slurm-1549075.out  uniquemappedelements.txt
RMHex02_bcassoc_allreads_barcodes_per_candidate-no_jackpots.feather             all_fastqjoin_merged.sorted.sam                              samfilefacts.err            slurm-1549076.out
RMHex02_bcassoc_allreads_barcodes_per_candidate-no_repeats-no_jackpots.feather  bamtosam3.sh                                                 samfilefacts.out            slurm-1549077.out
RMHex02_bcassoc_allreads_barcodes_per_candidate-no_repeats.feather              count_fastq.txt                                              samfilefacts_newdup.sh      slurm-1549079.out
(base) [rhauser@login02 MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing]$ nano 1_count_bc_nolab_new_try3.sh 
(base) [rhauser@login02 MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing]$ nano 3_PE_merge2.sh 
(base) [rhauser@login02 MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing]$ cd ..
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd MPRAflowBCassoc
(base) [rhauser@login02 MPRAflowBCassoc]$ ls
mpraflow.error  mpraflow.out  mpraflow_commentoutslurm.error  mpraflow_commentoutslurm.out  output  rmhex02bcassocmpraflow.sh  rmhex02bcassocmpraflow_slurmcommentout.sh  work
(base) [rhauser@login02 MPRAflowBCassoc]$ cd output/
(base) [rhauser@login02 output]$ ls
assoc_basic_RMHex02
(base) [rhauser@login02 output]$ cd assoc_basic_RMHex02/
(base) [rhauser@login02 assoc_basic_RMHex02]$ ls
count_fastq.txt  design_rmIllegalChars.fa  label_rmIllegalChars.txt
(base) [rhauser@login02 assoc_basic_RMHex02]$ cd ..
(base) [rhauser@login02 output]$ ls
assoc_basic_RMHex02
(base) [rhauser@login02 output]$ cd ..
(base) [rhauser@login02 MPRAflowBCassoc]$ ls
mpraflow.error  mpraflow.out  mpraflow_commentoutslurm.error  mpraflow_commentoutslurm.out  output  rmhex02bcassocmpraflow.sh  rmhex02bcassocmpraflow_slurmcommentout.sh  work
(base) [rhauser@login02 MPRAflowBCassoc]$ nano rmhex02bcassocmpraflow.sh
(base) [rhauser@login02 MPRAflowBCassoc]$ cd ..
(base) [rhauser@login02 rmhex02]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv          average_bedgraph4020698.out    finalrmhex02designfile.fa.ann                               onesidedneurons3repswithcontrolsneglog10pval.sorted.bed
CRISPRi                                                    average_bedgraph4020699.out    finalrmhex02designfile.fa.bwt                               onesidedneurons3repswithcontrolsneglog10pval.sorted.mean.bed
MPRAflowBCassoc                                            average_bedgraph4020700.out    finalrmhex02designfile.fa.dict                              output_Jane1
MPRAflowBCassoc_piecebypiece_notstrandspecific             average_bedgraph4020705.out    finalrmhex02designfile.fa.fai                               rmhex02_hekssignal.bed
MPRAflowBCassoc_piecebypiece_notstrandspecific_nodemuxing  barcode_association            finalrmhex02designfile.fa.pac                               rmhex02_hekssignal.col5bed.bed
MPRAflowCount                                              barcode_association_Jane       finalrmhex02designfile.fa.sa                                rmhex02_hekssignal.partition.bed
MPRAflowCount_allreps                                      design_initial                 finalrmhex02designfile_justscram.fa                         rmhex02_hekssignal.sorted.bed
MPRAflowCount_reps123                                      designaddon                    onesidedneurons3repswithcontrolsneglog10pval.bed            rmhex02_hekssignal.sorted.mean.bed
MPRAflow_piecebypiece                                      finalrmhex02designfile.fa      onesidedneurons3repswithcontrolsneglog10pval.col5bed.bed    seqcenter_rmhex02rep3_4_rmhex01rep3_20240405
allfastqs                                                  finalrmhex02designfile.fa.amb  onesidedneurons3repswithcontrolsneglog10pval.partition.bed  testEmpirical_alpha_Neurons_results_labeled_prefilt_allcontrols_oneside_alllables.tsv
(base) [rhauser@login02 rmhex02]$ cd ..
(base) [rhauser@login02 files]$ cd rmhex03
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd bcassoc/
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd bcassoc219
(base) [rhauser@login02 bcassoc219]$ ls
MPRAflow_assoc_219_7293347.out                    MPRAflow_assoc_219_called_shorttime2_mapq0cigar_7588835.out      MPRAflow_assoc_219ori_7294298.out                       outputcalled                  workingcalled2
MPRAflow_assoc_219_7294196.out                    MPRAflow_assoc_219_called_shorttime2_mapq0cigar_7592701.out      bcassoc219_edits2.sh                                    outputcalled2                 workingcalled2_mincov1
MPRAflow_assoc_219_7294264.out                    MPRAflow_assoc_219_called_shorttime2_mapq0cigar_7600330.out      bcassoc219_edits3.sh                                    outputcalled2_mincov1         workingcalled2_nocigarfilt
MPRAflow_assoc_219_called_7395166.out             MPRAflow_assoc_219_called_shorttime2_mincov1_7583170.out         bcassoc219_onlycalledreads.sh                           outputcalled2_nocigarfilt     workingcalled2_prepairedreads
MPRAflow_assoc_219_called_7439711.out             MPRAflow_assoc_219_called_shorttime2_prepairedreads_7619296.out  bcassoc219_onlycalledreads_mincov1.sh                   outputcalled2_prepairedreads  workingcalled_mapq0cigar
MPRAflow_assoc_219_called_shorttime2_7507659.out  MPRAflow_assoc_219edits2_7294336.out                             bcassoc219_onlycalledreadsshorttime2.sh                 outputcalled_mapq0cigar       workingchunk
MPRAflow_assoc_219_called_shorttime2_7507663.out  MPRAflow_assoc_219edits2_7316317.out                             bcassoc219_onlycalledreadsshorttime2_mapq0cigar.sh      outputchunk                   workingoriginalnf
MPRAflow_assoc_219_called_shorttime2_7557508.out  MPRAflow_assoc_219edits2_7316319.out                             bcassoc219_onlycalledreadsshorttime2_nocigarfilt.sh     outputoriginalnf
MPRAflow_assoc_219_called_shorttime2_7557512.out  MPRAflow_assoc_219edits2_7319941.out                             bcassoc219_onlycalledreadsshorttime2_prepairedreads.sh  working
MPRAflow_assoc_219_called_shorttime2_7580620.out  MPRAflow_assoc_219edits3_7320159.out                             output                                                  workingcalled
(base) [rhauser@login02 bcassoc219]$ nano bcassoc219_onlycalledreadsshorttime2_mapq0cigar.sh  
(base) [rhauser@login02 bcassoc219]$ cd ..
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd bcassoc218
(base) [rhauser@login02 bcassoc218]$ ls
MPRAflow_assoc_218_7293354.out  MPRAflow_assoc_218_called_shorttime2_7566671.out  bcassoc218.sh                            bcassoc218_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_218_7294464.out  MPRAflow_assoc_218_called_shorttime2_7614909.out  bcassoc218_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc218]$ nano bcassoc218_onlycalledreadsshorttime2_mapq0cigar.sh 
(base) [rhauser@login02 bcassoc218]$ cd ..
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd bcassoc216
(base) [rhauser@login02 bcassoc216]$ ls
MPRAflow_assoc_216_7293373.out  MPRAflow_assoc_216_called_shorttime2_7566675.out       bcassoc216.sh                            bcassoc216_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working
MPRAflow_assoc_216_7294472.out  MPRAflow_assoc_216_called_shorttime2_mapq_7614932.out  bcassoc216_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalledmapq
(base) [rhauser@login02 bcassoc216]$ nano bcassoc216_onlycalledreadsshorttime2_mapq0cigar.sh 
(base) [rhauser@login02 bcassoc216]$ cd ..
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd bcassoc214/
(base) [rhauser@login02 bcassoc214]$ ls
MPRAflow_assoc_214_7293429.out  MPRAflow_assoc_214_called_shorttime2_7566691.out  bcassoc214.sh                            bcassoc214_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_214_7294480.out  MPRAflow_assoc_214_called_shorttime2_7614930.out  bcassoc214_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc214]$ nano bcassoc214_onlycalledreadsshorttime2_mapq0cigar.sh 
(base) [rhauser@login02 bcassoc214]$ cd ..
(base) [rhauser@login02 bcassoc]$ cd bcassoc212
(base) [rhauser@login02 bcassoc212]$ ls
MPRAflow_assoc_212_7293461.out  MPRAflow_assoc_212_called_shorttime2_7566704.out      bcassoc212.sh                            bcassoc212_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_212_7294496.out  MPRAflow_assoc_212_called_shorttime2mapq_7614980.out  bcassoc212_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc212]$ nano bcassoc212_onlycalledreadsshorttime2_mapq0cigar.sh 
(base) [rhauser@login02 bcassoc212]$ cd ..
(base) [rhauser@login02 bcassoc]$ cd bcassoc210
(base) [rhauser@login02 bcassoc210]$ ls
MPRAflow_assoc_210_7293509.out  MPRAflow_assoc_210_called_shorttime2_7566722.out            bcassoc210.sh                            bcassoc210_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_210_7294512.out  MPRAflow_assoc_210_called_shorttime2_mapqcigar_7615014.out  bcassoc210_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc210]$ nano bcassoc210_onlycalledreadsshorttime2_mapq0cigar.sh 
(base) [rhauser@login02 bcassoc210]$ cd ..
(base) [rhauser@login02 bcassoc]$ cd bcassoc208/
(base) [rhauser@login02 bcassoc208]$ ls
MPRAflow_assoc_208_7293550.out  MPRAflow_assoc_208_called_shorttime2_7566785.out      bcassoc208.sh                            bcassoc208_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_208_7294529.out  MPRAflow_assoc_208_called_shorttime2mapq_7615054.out  bcassoc208_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc208]$ nano bcassoc208_onlycalledreadsshorttime2_mapq0cigar.sh
(base) [rhauser@login02 bcassoc208]$ cd ..
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd ..
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd design/
(base) [rhauser@login02 design]$ ls
ADSP_chr17_variation_mpra_overlap_9-5-24.list.tsv  fiveallelesmpradesign.tsv  mapt_vcf.tsv               sevenallelesmpradesign.tsv  variant_reference_key.tsv  varmpradesign_216bp.fasta      varmpradesign_219bp.fasta.pac       varmpradesign_219bp_refs.fasta.pac
combinedmpradesignfile.tsv                         fourallelesmpradesign.tsv  mapt_vcf.tsv.gz            sizecontrolstest.tsv        variantsfromjared          varmpradesign_218bp.fasta      varmpradesign_219bp.fasta.sa        varmpradesign_219bp_refs.fasta.sa
cutvcf.sh                                          indels_vcf.tsv             mapt_vcf.vcf               snvs_vcf.tsv                varmpradesign_208bp.fasta  varmpradesign_219bp.fasta      varmpradesign_219bp_refs.fasta      varmpradesign_all.fasta
cutvsf.out                                         indels_vcf.tsv.gz          motifbreakersetupsnps.bed  snvs_vcf.tsv.gz             varmpradesign_210bp.fasta  varmpradesign_219bp.fasta.amb  varmpradesign_219bp_refs.fasta.amb
design_allsizes_rmIllegalChars.fa                  indels_vcf.vcf             motifbreakersetupsnps.tsv  snvs_vcf.vcf                varmpradesign_212bp.fasta  varmpradesign_219bp.fasta.ann  varmpradesign_219bp_refs.fasta.ann
designscript                                       labelfile_rmhex03.tsv      newfinalvarmpra2.fasta     twoallelesmpradesign.tsv    varmpradesign_214bp.fasta  varmpradesign_219bp.fasta.bwt  varmpradesign_219bp_refs.fasta.bwt
(base) [rhauser@login02 design]$ cd ..
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd bcassoc/
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd bcassoc210
(base) [rhauser@login02 bcassoc210]$ ls
MPRAflow_assoc_210_7293509.out  MPRAflow_assoc_210_called_shorttime2_7566722.out            bcassoc210.sh                            bcassoc210_onlycalledreadsshorttime2_mapq0cigar.sh  outputcalled2     working         workingcalledmapq
MPRAflow_assoc_210_7294512.out  MPRAflow_assoc_210_called_shorttime2_mapqcigar_7615014.out  bcassoc210_onlycalledreadsshorttime2.sh  output                                              outputcalledmapq  workingcalled2
(base) [rhauser@login02 bcassoc210]$ cd outputcalledmapq/
(base) [rhauser@login02 outputcalledmapq]$ ls
assoc_basic_210_called_shorttime2_mapq
(base) [rhauser@login02 outputcalledmapq]$ cd assoc_basic_210_called_shorttime2_mapq/
(base) [rhauser@login02 assoc_basic_210_called_shorttime2_mapq]$ ls
assoc_basic_210_called_shorttime2_mapq_barcode_counts.pickle                                  assoc_basic_210_called_shorttime2_mapq_filtered_coords_to_barcodes.tsv  design_rmIllegalChars.fa    original_count_summary.txt
assoc_basic_210_called_shorttime2_mapq_barcodes_per_candidate-no_repeats-no_jackpots.feather  assoc_basic_210_called_shorttime2_mapq_filtered_counts.png              dicttotsv.sh                slurm-7621093.out
assoc_basic_210_called_shorttime2_mapq_coords_to_barcodes.pickle                              assoc_basic_210_called_shorttime2_mapq_original_counts.png              filtered_count_summary.txt
assoc_basic_210_called_shorttime2_mapq_filtered_coords_to_barcodes.pickle                     count_fastq.txt                                                         label_rmIllegalChars.txt
(base) [rhauser@login02 assoc_basic_210_called_shorttime2_mapq]$ nano dicttotsv.sh 
(base) [rhauser@login02 assoc_basic_210_called_shorttime2_mapq]$ nano dicttotsv.sh 
(base) [rhauser@login02 assoc_basic_210_called_shorttime2_mapq]$ nano /cluster/home/rhauser/scripts/mpraflowdict_to_tsv.py
(base) [rhauser@login02 assoc_basic_210_called_shorttime2_mapq]$ cd ..
(base) [rhauser@login02 outputcalledmapq]$ ls
assoc_basic_210_called_shorttime2_mapq
(base) [rhauser@login02 outputcalledmapq]$ cd ..
(base) [rhauser@login02 bcassoc210]$ cd ..
(base) [rhauser@login02 bcassoc]$ ls
barcodeAssociationSeqDemux_noLaneSplitting.sh             bcassoc208  bcassoc212  bcassoc216  bcassoc219      bcassoc_reseq  demux.out            seq
barcodeAssociationSeqDemux_noLaneSplitting_reseq_trim.sh  bcassoc210  bcassoc214  bcassoc218  bcassoc219refs  bwa219test     demux_reseqtrim.out
(base) [rhauser@login02 bcassoc]$ cd ..
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd MPRAsnakeflowcount/
(base) [rhauser@login02 MPRAsnakeflowcount]$ ls
MPRAsnakeflow_count_7715787.out  MPRAsnakeflow_count_7716197.out  MPRAsnakeflow_count_7717566.out  MPRAsnakeflow_count_7717689.out  MPRAsnakeflow_count_7986016.out  MPRAsnakeflow_count_7986050.out  count_justneurons.yaml                   mpranalyze
MPRAsnakeflow_count_7715788.out  MPRAsnakeflow_count_7716198.out  MPRAsnakeflow_count_7717571.out  MPRAsnakeflow_count_7717704.out  MPRAsnakeflow_count_7986024.out  bcalm                            experiment_file_rmhex03.csv              mprasnakeflowcall.sh
MPRAsnakeflow_count_7715789.out  MPRAsnakeflow_count_7716199.out  MPRAsnakeflow_count_7717610.out  MPRAsnakeflow_count_7717753.out  MPRAsnakeflow_count_7986030.out  count.yaml                       experiment_file_rmhex03_justheks.csv     mprasnakeflowcalltest.sh
MPRAsnakeflow_count_7715790.out  MPRAsnakeflow_count_7716200.out  MPRAsnakeflow_count_7717624.out  MPRAsnakeflow_count_7718361.out  MPRAsnakeflow_count_7986038.out  count_justheks.yaml              experiment_file_rmhex03_justneurons.csv  mprasnakeflowcalltest2.sh
MPRAsnakeflow_count_7716195.out  MPRAsnakeflow_count_7716203.out  MPRAsnakeflow_count_7717674.out  MPRAsnakeflow_count_7986011.out  MPRAsnakeflow_count_7986047.out  count_justheks2.yaml             experiment_file_rmhex03test.csv          results
(base) [rhauser@login02 MPRAsnakeflowcount]$ nano mprasnakeflowcalltest2.sh
(base) [rhauser@login02 MPRAsnakeflowcount]$ ls
MPRAsnakeflow_count_7715787.out  MPRAsnakeflow_count_7716197.out  MPRAsnakeflow_count_7717566.out  MPRAsnakeflow_count_7717689.out  MPRAsnakeflow_count_7986016.out  MPRAsnakeflow_count_7986050.out  count_justneurons.yaml                   mpranalyze
MPRAsnakeflow_count_7715788.out  MPRAsnakeflow_count_7716198.out  MPRAsnakeflow_count_7717571.out  MPRAsnakeflow_count_7717704.out  MPRAsnakeflow_count_7986024.out  bcalm                            experiment_file_rmhex03.csv              mprasnakeflowcall.sh
MPRAsnakeflow_count_7715789.out  MPRAsnakeflow_count_7716199.out  MPRAsnakeflow_count_7717610.out  MPRAsnakeflow_count_7717753.out  MPRAsnakeflow_count_7986030.out  count.yaml                       experiment_file_rmhex03_justheks.csv     mprasnakeflowcalltest.sh
MPRAsnakeflow_count_7715790.out  MPRAsnakeflow_count_7716200.out  MPRAsnakeflow_count_7717624.out  MPRAsnakeflow_count_7718361.out  MPRAsnakeflow_count_7986038.out  count_justheks.yaml              experiment_file_rmhex03_justneurons.csv  mprasnakeflowcalltest2.sh
MPRAsnakeflow_count_7716195.out  MPRAsnakeflow_count_7716203.out  MPRAsnakeflow_count_7717674.out  MPRAsnakeflow_count_7986011.out  MPRAsnakeflow_count_7986047.out  count_justheks2.yaml             experiment_file_rmhex03test.csv          results
(base) [rhauser@login02 MPRAsnakeflowcount]$ nano count.yaml
(base) [rhauser@login02 MPRAsnakeflowcount]$ cd ~
(base) [rhauser@login02 ~]$ cd MPRAsnakeflow/
(base) [rhauser@login02 MPRAsnakeflow]$ ls
CHANGELOG.md  LICENSE    config  logs            profiles   results               rmhex01snakemake_2.log         rmhex03snakemake_justheks2.log    rmhex03snakemake_main.log   rmhex03snakemake_main3.log  workflow
Dockerfile    README.md  docs    maptforbri.log  resources  rmhex01snakemake.log  rmhex03snakemake_justheks.log  rmhex03snakemake_justneurons.log  rmhex03snakemake_main2.log  version.txt
(base) [rhauser@login02 MPRAsnakeflow]$ cd profiles/
(base) [rhauser@login02 profiles]$ ls
config_maptforbri  default  default_2  default_3  default_rmh  results  rmhex01_redo  rmhex01_redo2  rmhex01_redo3  rmhex03_neurons  rmhex03snakemake_justheks.log
(base) [rhauser@login02 profiles]$ cd ..
(base) [rhauser@login02 MPRAsnakeflow]$ ls
CHANGELOG.md  LICENSE    config  logs            profiles   results               rmhex01snakemake_2.log         rmhex03snakemake_justheks2.log    rmhex03snakemake_main.log   rmhex03snakemake_main3.log  workflow
Dockerfile    README.md  docs    maptforbri.log  resources  rmhex01snakemake.log  rmhex03snakemake_justheks.log  rmhex03snakemake_justneurons.log  rmhex03snakemake_main2.log  version.txt
(base) [rhauser@login02 MPRAsnakeflow]$ cd ~/files/rmhex03
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd MPRAsnakeflowcount/
(base) [rhauser@login02 MPRAsnakeflowcount]$ ls
MPRAsnakeflow_count_7715787.out  MPRAsnakeflow_count_7716197.out  MPRAsnakeflow_count_7717566.out  MPRAsnakeflow_count_7717689.out  MPRAsnakeflow_count_7986016.out  MPRAsnakeflow_count_7986050.out  count_justneurons.yaml                   mpranalyze
MPRAsnakeflow_count_7715788.out  MPRAsnakeflow_count_7716198.out  MPRAsnakeflow_count_7717571.out  MPRAsnakeflow_count_7717704.out  MPRAsnakeflow_count_7986024.out  bcalm                            experiment_file_rmhex03.csv              mprasnakeflowcall.sh
MPRAsnakeflow_count_7715789.out  MPRAsnakeflow_count_7716199.out  MPRAsnakeflow_count_7717610.out  MPRAsnakeflow_count_7717753.out  MPRAsnakeflow_count_7986030.out  count.yaml                       experiment_file_rmhex03_justheks.csv     mprasnakeflowcalltest.sh
MPRAsnakeflow_count_7715790.out  MPRAsnakeflow_count_7716200.out  MPRAsnakeflow_count_7717624.out  MPRAsnakeflow_count_7718361.out  MPRAsnakeflow_count_7986038.out  count_justheks.yaml              experiment_file_rmhex03_justneurons.csv  mprasnakeflowcalltest2.sh
MPRAsnakeflow_count_7716195.out  MPRAsnakeflow_count_7716203.out  MPRAsnakeflow_count_7717674.out  MPRAsnakeflow_count_7986011.out  MPRAsnakeflow_count_7986047.out  count_justheks2.yaml             experiment_file_rmhex03test.csv          results
(base) [rhauser@login02 MPRAsnakeflowcount]$ cd bcalm/
(base) [rhauser@login02 bcalm]$ ls
complexheatmap.pdf           indelcorr_logFCcutoff_10bcthres.pdf  maptbcalm_logfc_30bases.pdf                  maptvol_neurons_10bcthres.pdf              snvsvol_heks_filt_10bcthres_logFCcutoff.pdf
heks                         indelvol_heks_filt_10bcthres.pdf     maptbcalm_logfc_50basesdot.pdf               maptvol_neurons_10bcthres_logFCcutoff.pdf  snvsvol_heks_nofilt.pdf
indelcorr.pdf                indelvol_heks_nofilt.pdf             maptheatmap_bcalm_10bcthres.pdf              neurons                                    snvsvol_neurons_filt_10bcthres.pdf
indelcorr_filt10bcthres.pdf  indelvol_neurons_filt_10bcthres.pdf  maptheatmap_bcalm_10bcthres_logFCcutoff.pdf  snvcorr_logFCcutoff_10bcthres.pdf          snvsvol_neurons_filt_10bcthres_logFCcutoff.pdf
indelcorr_logFCcutoff.pdf    indelvol_neurons_nofilt.pdf          maptvol_neurons.pdf                          snvsvol_heks_filt_10bcthres.pdf            snvsvol_neurons_nofilt.pdf
(base) [rhauser@login02 bcalm]$ cd heks/
(base) [rhauser@login02 heks]$ ls
bcalm_8038629.out  bcalm_8038721.out  bcalm_8038809.out  bcalm_8039460.out  bcalm_8047207.out              bcalm_indels_labeled_heks_10bcthres.tsv             bcalm_snvsandindels_labeled_heks.tsv            bcalmmapt_10bc.R    bcalmsnvsindels_10bc.R
bcalm_8038631.out  bcalm_8038747.out  bcalm_8038818.out  bcalm_8039464.out  bcalm_8052883.out              bcalm_indels_labeled_heks_10bcthres_fulllabels.tsv  bcalm_snvsandindels_labeled_heks_10bcthres.tsv  bcalmmapt_10bc.sh   bcalmsnvsindels_10bc.sh
bcalm_8038633.out  bcalm_8038754.out  bcalm_8038824.out  bcalm_8039466.out  bcalm_8055757.out              bcalm_mapt_heks.tsv                                 bcalmindel.R                                    bcalmsnvs.R
bcalm_8038672.out  bcalm_8038755.out  bcalm_8038831.out  bcalm_8039467.out  bcalm_8055758.out              bcalm_mapt_heks_10bcsthres.tsv                      bcalmindel.sh                                   bcalmsnvs.sh
bcalm_8038680.out  bcalm_8038762.out  bcalm_8038832.out  bcalm_8044702.out  bcalm_8055759.out              bcalm_mapt_heks_10bcsthres_labeled.tsv              bcalmindel_10bc.R                               bcalmsnvs_10bc.R
bcalm_8038707.out  bcalm_8038784.out  bcalm_8038884.out  bcalm_8044828.out  bcalm_8055760.out              bcalm_snvs_labeled_heks.tsv                         bcalmindel_10bc.sh                              bcalmsnvs_10bc.sh
bcalm_8038719.out  bcalm_8038797.out  bcalm_8038885.out  bcalm_8044829.out  bcalm_8082996.out              bcalm_snvs_labeled_heks_10bcthres.tsv               bcalmmapt.R                                     bcalmsnvsindels.R
bcalm_8038720.out  bcalm_8038801.out  bcalm_8038886.out  bcalm_8047206.out  bcalm_indels_labeled_heks.tsv  bcalm_snvs_labeled_heks_10bcthres_fulllabels.tsv    bcalmmapt.sh                                    bcalmsnvsindels.sh
(base) [rhauser@login02 heks]$ nano bcalmmapt.R
(base) [rhauser@login02 heks]$ cd ..
(base) [rhauser@login02 bcalm]$ ls
complexheatmap.pdf           indelcorr_logFCcutoff_10bcthres.pdf  maptbcalm_logfc_30bases.pdf                  maptvol_neurons_10bcthres.pdf              snvsvol_heks_filt_10bcthres_logFCcutoff.pdf
heks                         indelvol_heks_filt_10bcthres.pdf     maptbcalm_logfc_50basesdot.pdf               maptvol_neurons_10bcthres_logFCcutoff.pdf  snvsvol_heks_nofilt.pdf
indelcorr.pdf                indelvol_heks_nofilt.pdf             maptheatmap_bcalm_10bcthres.pdf              neurons                                    snvsvol_neurons_filt_10bcthres.pdf
indelcorr_filt10bcthres.pdf  indelvol_neurons_filt_10bcthres.pdf  maptheatmap_bcalm_10bcthres_logFCcutoff.pdf  snvcorr_logFCcutoff_10bcthres.pdf          snvsvol_neurons_filt_10bcthres_logFCcutoff.pdf
indelcorr_logFCcutoff.pdf    indelvol_neurons_nofilt.pdf          maptvol_neurons.pdf                          snvsvol_heks_filt_10bcthres.pdf            snvsvol_neurons_nofilt.pdf
(base) [rhauser@login02 bcalm]$ cd ..
(base) [rhauser@login02 MPRAsnakeflowcount]$ ls
MPRAsnakeflow_count_7715787.out  MPRAsnakeflow_count_7716197.out  MPRAsnakeflow_count_7717566.out  MPRAsnakeflow_count_7717689.out  MPRAsnakeflow_count_7986016.out  MPRAsnakeflow_count_7986050.out  count_justneurons.yaml                   mpranalyze
MPRAsnakeflow_count_7715788.out  MPRAsnakeflow_count_7716198.out  MPRAsnakeflow_count_7717571.out  MPRAsnakeflow_count_7717704.out  MPRAsnakeflow_count_7986024.out  bcalm                            experiment_file_rmhex03.csv              mprasnakeflowcall.sh
MPRAsnakeflow_count_7715789.out  MPRAsnakeflow_count_7716199.out  MPRAsnakeflow_count_7717610.out  MPRAsnakeflow_count_7717753.out  MPRAsnakeflow_count_7986030.out  count.yaml                       experiment_file_rmhex03_justheks.csv     mprasnakeflowcalltest.sh
MPRAsnakeflow_count_7715790.out  MPRAsnakeflow_count_7716200.out  MPRAsnakeflow_count_7717624.out  MPRAsnakeflow_count_7718361.out  MPRAsnakeflow_count_7986038.out  count_justheks.yaml              experiment_file_rmhex03_justneurons.csv  mprasnakeflowcalltest2.sh
MPRAsnakeflow_count_7716195.out  MPRAsnakeflow_count_7716203.out  MPRAsnakeflow_count_7717674.out  MPRAsnakeflow_count_7986011.out  MPRAsnakeflow_count_7986047.out  count_justheks2.yaml             experiment_file_rmhex03test.csv          results
(base) [rhauser@login02 MPRAsnakeflowcount]$ cd ..
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd ..
(base) [rhauser@login02 files]$ ls
andersonrogers2023supplementaldata  encode  fastqlist.txt  from_henry  from_nick  from_others  graphs  hic_rizzardi  practice  rmhex01  rmhex02  rmhex03  rogers2024supplementaldata  slurm-8238174.out  slurm-8238956.out  temp
(base) [rhauser@login02 files]$ cd ..
(base) [rhauser@login02 ~]$ cd MPRAsnakeflow/
(base) [rhauser@login02 MPRAsnakeflow]$ ls
CHANGELOG.md  LICENSE    config  logs            profiles   results               rmhex01snakemake_2.log         rmhex03snakemake_justheks2.log    rmhex03snakemake_main.log   rmhex03snakemake_main3.log  workflow
Dockerfile    README.md  docs    maptforbri.log  resources  rmhex01snakemake.log  rmhex03snakemake_justheks.log  rmhex03snakemake_justneurons.log  rmhex03snakemake_main2.log  version.txt
(base) [rhauser@login02 MPRAsnakeflow]$ nano rmhex03snakemake_main3.log
(base) [rhauser@login02 MPRAsnakeflow]$ cd results/
(base) [rhauser@login02 results]$ ls
assignment  experiments  logs
(base) [rhauser@login02 results]$ cd experiments/
(base) [rhauser@login02 experiments]$ ls
rmhex01Count  rmhex01Count_2  rmhex03Count  rmhex03Count_justheks2  rmhex03Count_justneurons
(base) [rhauser@login02 experiments]$ cd rmhex03Count
(base) [rhauser@login02 rmhex03Count]$ ls
assigned_counts                                            reporter_experiment.barcode.HEKs.fromFile.default.all.tsv.gz                        reporter_experiment.oligo.HEKs.fromFile.default.all.tsv.gz
assignment                                                 reporter_experiment.barcode.HEKs.fromFile.default.min_oligo_threshold_10.tsv        reporter_experiment.oligo.HEKs.fromFile.default.min_oligo_threshold_10.tsv.gz
qc_metrics.HEKs.fromFile.default.json                      reporter_experiment.barcode.HEKs.fromFile.default.min_oligo_threshold_10.tsv.gz     reporter_experiment.oligo.Neurons.fromFile.default.all.tsv.gz
qc_metrics.Neurons.fromFile.default.json                   reporter_experiment.barcode.Neurons.fromFile.default.all.tsv                        reporter_experiment.oligo.Neurons.fromFile.default.min_oligo_threshold_10.tsv.gz
qc_report.HEKs.fromFile.default.html                       reporter_experiment.barcode.Neurons.fromFile.default.all.tsv.gz                     statistic
qc_report.Neurons.fromFile.default.html                    reporter_experiment.barcode.Neurons.fromFile.default.min_oligo_threshold_10.tsv
reporter_experiment.barcode.HEKs.fromFile.default.all.tsv  reporter_experiment.barcode.Neurons.fromFile.default.min_oligo_threshold_10.tsv.gz
(base) [rhauser@login02 rmhex03Count]$ cd ..
(base) [rhauser@login02 experiments]$ cd ..
(base) [rhauser@login02 results]$ ls
assignment  experiments  logs
(base) [rhauser@login02 results]$ cd ..
(base) [rhauser@login02 MPRAsnakeflow]$ ls
CHANGELOG.md  LICENSE    config  logs            profiles   results               rmhex01snakemake_2.log         rmhex03snakemake_justheks2.log    rmhex03snakemake_main.log   rmhex03snakemake_main3.log  workflow
Dockerfile    README.md  docs    maptforbri.log  resources  rmhex01snakemake.log  rmhex03snakemake_justheks.log  rmhex03snakemake_justneurons.log  rmhex03snakemake_main2.log  version.txt
(base) [rhauser@login02 MPRAsnakeflow]$ cd profiles/
(base) [rhauser@login02 profiles]$ ls
config_maptforbri  default  default_2  default_3  default_rmh  results  rmhex01_redo  rmhex01_redo2  rmhex01_redo3  rmhex03_neurons  rmhex03snakemake_justheks.log
(base) [rhauser@login02 profiles]$ cd default_3
(base) [rhauser@login02 default_3]$ ls
config.yaml
(base) [rhauser@login02 default_3]$ nano config.yaml 
(base) [rhauser@login02 default_3]$ cd ..
(base) [rhauser@login02 profiles]$ cd default_2
(base) [rhauser@login02 default_2]$ ls
config.yaml
(base) [rhauser@login02 default_2]$ nano config.yaml 
(base) [rhauser@login02 default_2]$ cd ..
(base) [rhauser@login02 profiles]$ cd default
(base) [rhauser@login02 default]$ ls
config.yaml
(base) [rhauser@login02 default]$ nano config.yaml 
(base) [rhauser@login02 default]$ cd ..
(base) [rhauser@login02 profiles]$ ls
config_maptforbri  default  default_2  default_3  default_rmh  results  rmhex01_redo  rmhex01_redo2  rmhex01_redo3  rmhex03_neurons  rmhex03snakemake_justheks.log
(base) [rhauser@login02 profiles]$ cd ..
(base) [rhauser@login02 MPRAsnakeflow]$ ls
CHANGELOG.md  LICENSE    config  logs            profiles   results               rmhex01snakemake_2.log         rmhex03snakemake_justheks2.log    rmhex03snakemake_main.log   rmhex03snakemake_main3.log  workflow
Dockerfile    README.md  docs    maptforbri.log  resources  rmhex01snakemake.log  rmhex03snakemake_justheks.log  rmhex03snakemake_justneurons.log  rmhex03snakemake_main2.log  version.txt
(base) [rhauser@login02 MPRAsnakeflow]$ cd ~
(base) [rhauser@login02 ~]$ cd files/rmhex03
(base) [rhauser@login02 rmhex03]$ ls
AIfromHenry  MPRAflowCount  MPRAflowCount_newSeq  MPRAsnakeflowcount  bcassoc  combinedFastqs  design  forGEO  meme  motifbreakr  newOutputSequencing  outputSequencing  spliceAI  variantinfo
(base) [rhauser@login02 rmhex03]$ cd MPRAsnakeflowcount/
(base) [rhauser@login02 MPRAsnakeflowcount]$ ls
MPRAsnakeflow_count_7715787.out  MPRAsnakeflow_count_7716197.out  MPRAsnakeflow_count_7717566.out  MPRAsnakeflow_count_7717689.out  MPRAsnakeflow_count_7986016.out  MPRAsnakeflow_count_7986050.out  count_justneurons.yaml                   mpranalyze
MPRAsnakeflow_count_7715788.out  MPRAsnakeflow_count_7716198.out  MPRAsnakeflow_count_7717571.out  MPRAsnakeflow_count_7717704.out  MPRAsnakeflow_count_7986024.out  bcalm                            experiment_file_rmhex03.csv              mprasnakeflowcall.sh
MPRAsnakeflow_count_7715789.out  MPRAsnakeflow_count_7716199.out  MPRAsnakeflow_count_7717610.out  MPRAsnakeflow_count_7717753.out  MPRAsnakeflow_count_7986030.out  count.yaml                       experiment_file_rmhex03_justheks.csv     mprasnakeflowcalltest.sh
MPRAsnakeflow_count_7715790.out  MPRAsnakeflow_count_7716200.out  MPRAsnakeflow_count_7717624.out  MPRAsnakeflow_count_7718361.out  MPRAsnakeflow_count_7986038.out  count_justheks.yaml              experiment_file_rmhex03_justneurons.csv  mprasnakeflowcalltest2.sh
MPRAsnakeflow_count_7716195.out  MPRAsnakeflow_count_7716203.out  MPRAsnakeflow_count_7717674.out  MPRAsnakeflow_count_7986011.out  MPRAsnakeflow_count_7986047.out  count_justheks2.yaml             experiment_file_rmhex03test.csv          results
(base) [rhauser@login02 MPRAsnakeflowcount]$ cd bcalm/
(base) [rhauser@login02 bcalm]$ ls
complexheatmap.pdf           indelcorr_logFCcutoff_10bcthres.pdf  maptbcalm_logfc_30bases.pdf                  maptvol_neurons_10bcthres.pdf              snvsvol_heks_filt_10bcthres_logFCcutoff.pdf
heks                         indelvol_heks_filt_10bcthres.pdf     maptbcalm_logfc_50basesdot.pdf               maptvol_neurons_10bcthres_logFCcutoff.pdf  snvsvol_heks_nofilt.pdf
indelcorr.pdf                indelvol_heks_nofilt.pdf             maptheatmap_bcalm_10bcthres.pdf              neurons                                    snvsvol_neurons_filt_10bcthres.pdf
indelcorr_filt10bcthres.pdf  indelvol_neurons_filt_10bcthres.pdf  maptheatmap_bcalm_10bcthres_logFCcutoff.pdf  snvcorr_logFCcutoff_10bcthres.pdf          snvsvol_neurons_filt_10bcthres_logFCcutoff.pdf
indelcorr_logFCcutoff.pdf    indelvol_neurons_nofilt.pdf          maptvol_neurons.pdf                          snvsvol_heks_filt_10bcthres.pdf            snvsvol_neurons_nofilt.pdf
(base) [rhauser@login02 bcalm]$ cd neurons/
(base) [rhauser@login02 neurons]$ ls
bcalm_8038890.out  bcalm_indels_labeled_neurons.tsv                       bcalm_snvs_labeled_neurons_10bcthres_fulllabels.tsv  bcalmmapt.sh       bcalmsnvsindels.sh                       bricrispr_45848303_45849886_diamond.pdf  crispr_region9_pgdiamond.pdf
bcalm_8038891.out  bcalm_indels_labeled_neurons_10bcthres_fulllabels.tsv  bcalm_snvsandindels_labeled_neurons.tsv              bcalmmapt_10bc.R   bcalmsnvsindels_10bc.sh                  bricrispr_45873488_45874156_diamond.pdf  mapt_indels_diamonds.pdf
bcalm_8038892.out  bcalm_indels_labeled_neuronss_10bcthres.tsv            bcalm_snvsandindels_labeled_neurons_10bcthres.tsv    bcalmmapt_10bc.sh  bcalmsnvsindels_10bcthres.R              bricrispr_45940685_45941352_diamond.pdf  mpranalyze.R
bcalm_8038893.out  bcalm_mapt_neurons.tsv                                 bcalmindel.R                                         bcalmsnvs.R        bricrispr_45219148_45219936_diamond.pdf  bricrispr_45942082_45942750_diamond.pdf  mpranalyze.sh
bcalm_8038894.out  bcalm_mapt_neurons_10bcthres.tsv                       bcalmindel.sh                                        bcalmsnvs.sh       bricrispr_45240986_45242338_diamond.pdf  bricrispr_45971424_45972091_diamond.pdf  newcombined.tsv
bcalm_8038895.out  bcalm_mapt_neurons_10bcthres_withlabels.tsv            bcalmindel_10bc.R                                    bcalmsnvs_10bc.R   bricrispr_45309210_45309877_diamond.pdf  crispr_region18_PGscale_diamond.pdf      newcombined_onlymaptcres.tsv
bcalm_8044830.out  bcalm_snvs_labeled_neurons.tsv                         bcalmindel_10bc.sh                                   bcalmsnvs_10bc.sh  bricrispr_45428787_45429857_diamond.pdf  crispr_region5_PGscale_diamond.pdf
bcalm_8044831.out  bcalm_snvs_labeled_neurons_10bcthres.tsv               bcalmmapt.R                                          bcalmsnvsindels.R  bricrispr_45431444_45432656_diamond.pdf  crispr_region8_pgdiamond.pdf
(base) [rhauser@login02 neurons]$ nano bcalmindel_10bc.R
(base) [rhauser@login02 neurons]$ nano bcalmsnvs_10bc.R
(base) [rhauser@login02 neurons]$ nano bcalmmapt_10bc.R

  GNU nano 5.6.1                                                                                                               bcalmmapt_10bc.R                                                                                                                          
key <- read_tsv("/cluster/home/rhauser/files/rmhex03/design/fiveallelesmpradesign.tsv")

#filter key for just mapt
maptkey <- key %>% filter(str_detect(id, regex("mapt", ignore_case = TRUE)))
maptkey <- maptkey %>% filter(!str_detect(id, regex("scram", ignore_case = TRUE)))
newkey <- maptkey %>%
  group_by(element) %>%
  summarize(
    ID   = dplyr::first(element),
    REF  = id[alleleA],
    ALT1 = id[alleleB],
    ALT2 = id[alleleC],
    ALT3 = id[alleleD],
    ALT4 = id[alleleE],
    .groups = "drop"
  )


create_var_df_satmut <- function(df, map_df) {
  # Check required columns
  required_cols <- c("ID", "REF", "ALT1", "ALT2", "ALT3", "ALT4")
  if (!all(required_cols %in% colnames(map_df))) {
    stop("map_df must contain columns: 'ID', 'REF', 'ALT1', 'ALT2', 'ALT3', 'ALT4'")
  }

  if (!"name" %in% colnames(df)) {
    stop("df must contain column 'name'")
  }

  # Reshape map_df to long format for ALTs
  alt_cols <- c("ALT1", "ALT2", "ALT3", "ALT4")
  map_alt_long <- map_df %>%
    pivot_longer(cols = all_of(alt_cols), names_to = "alt_type", values_to = "ALT")

  # Check for matching refs or alts
  if (!any(df$name %in% map_df$REF) & !any(df$name %in% map_alt_long$ALT)) {
    stop("No matches found between the 'name' column in 'df' and the 'REF'/ALT columns in 'map_df'.")
  }

  # Merge on REF
  df_ref <- merge(df, map_df[, c("ID", "REF")], by.x = "name", by.y = "REF", all.x = FALSE)
  df_ref$allele <- "ref"

  # Merge on ALT (long format)
  df_alt <- merge(df, map_alt_long[, c("ID", "ALT", "alt_type")], by.x = "name", by.y = "ALT", all.x = FALSE)
  df_alt$allele <- df_alt$alt_type
  df_alt$alt_type <- NULL

  # Remove ALT from ref to avoid confusion
# Remove ALT from ref to avoid confusion
  df_ref$REF <- NULL

  # Combine results
  df_combined <- dplyr::bind_rows(df_ref, df_alt)

  # Select and return final output
  var_df <- df_combined %>%
    dplyr::select(variant_id = ID, allele, Barcode, dplyr::matches("count"))

  return(var_df)
}

vardf <- create_var_df_satmut(maptcounts, newkey)

dna_var <- create_dna_df(vardf)
rna_var <- create_rna_df(vardf)

mympra <- MPRASet(DNA = dna_var, RNA = rna_var, eid = row.names(dna_var), barcode = NULL)

nreps <- 4
bcs <- ncol(dna_var) / nreps
design <- data.frame(intcpt = 1, alt1 = grepl("ALT1", colnames(mympra)), alt2 = grepl("ALT2", colnames(mympra)), alt3 = grepl("ALT3", colnames(mympra)), alt4 = grepl("ALT4", colnames(mympra)))
block_vector <- rep(1:nreps, each=bcs)
mpralm_fit_var <- mpralm(object = mympra, design = design, aggregate = "none", normalize = TRUE, model_type = "corr_groups", plot = FALSE, block = block_vector)

top_var_1 <- topTable(mpralm_fit_var, coef = 2, number = Inf)
top_var_2 <- topTable(mpralm_fit_var, coef = 3, number = Inf)
top_var_3 <- topTable(mpralm_fit_var, coef = 4, number = Inf)
top_var_4 <- topTable(mpralm_fit_var, coef = 5, number = Inf)

top_var_1 <- top_var_1 %>%
  rownames_to_column(var="element")
top_var_2 <- top_var_2 %>%
  rownames_to_column(var="element")
top_var_3 <- top_var_3 %>%
  rownames_to_column(var="element")
top_var_4 <- top_var_4 %>%
  rownames_to_column(var="element")

top_var_1$allele <- "B"
top_var_2$allele <- "C"
top_var_3$allele <- "D"
top_var_4$allele <- "E"
toptab_mapt <- rbind(top_var_1, top_var_2, top_var_3, top_var_4)

#add element names back in


write_tsv(toptab_mapt, "bcalm_mapt_neurons_10bcthres.tsv")
