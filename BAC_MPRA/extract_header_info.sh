#extract_header_info_forBecky_2022November03.sh


################################################################################
################################################################################
#This Script will extract information from your FASTQs.  It expects a paired FASTQ input
#where there is a repeating 4 line structure in each file.  It expects the 1st line to be
#structured as:
#SeqName [read identifier][Barcode]+[IndexRead]
#And the second line to be the sequence of the read.
#It will then output a table with the following information, separated by tabs
#ReadName, Barcode1, Index1, Barcode2, Index2
#
#This table can then be used later to easily 1) restrict to only reads with
#consistent Barcode and Index reads, and compare header names with a .sam file
#to identify where each read pair maps to.
#We will use this information to relate a given Barcode to a region of the genome.
#
#The script is called using the following command:
#bash /path/to/script/extract_header_info_forBecky_2022November03.sh /path/to/desired/outDir /path/to/fastq1.fastq /path/to/fastq2.fastq
#
#
#test use for my purposes:
#bash /cluster/home/bmoyers/for_becky/testing_header_script/extract_header_info_forBecky_2022November03.sh /cluster/home/bmoyers/for_becky/testing_header_script /cluster/home/bmoyers/for_becky/testing_header_script/head_HT7GJBGXL_s1_1_BC_X.fastq /cluster/home/bmoyers/for_becky/testing_header_script/head_HT7GJBGXL_s1_2_BC_X.fastq
################################################################################
################################################################################

################################################################################
################################################################################
#Begin Script
################################################################################
################################################################################


################################################################################
#The following reads in the arguments you provided.
################################################################################

outDir=${1}
fastq1=${2}
fastq2=${3}

#outDir="/cluster/home/bmoyers/for_becky/testing_header_script"
#fastq1="/cluster/home/bmoyers/for_becky/testing_header_script/head_HT7GJBGXL_s1_1_BC_X.fastq"
#fastq2="/cluster/home/bmoyers/for_becky/testing_header_script/head_HT7GJBGXL_s1_2_BC_X.fastq"

################################################################################
#The following couple of lines are fancy bits of script that will
#just pull out a basename from the provided fastq names.
#For each filename, it removes all of the name up to and including the last "/"
#and then removes everything after the last "." and beyond.
#So if the input is:
#"/cluster/home/rhauser/files/rmhex01/barcodeassocfastqs/HT7GJBGXL_s1_1_BC_X.fastq"
#Then the basename will become "HT7GJBGXL_s1_1_BC_X"
#It also creates the final baseName, which just shaves off the last "_" and everything
#behind it 3 times.  So this would leave, in the above example:
#HT7GJBGXL_s1
#There may be a more eloquent way to do this.  If we start getting weird
#output because the expected input filename format is different, we may have to rewrite.
################################################################################

baseName1="${fastq1##*/}"
baseName1="${baseName1%.*}"

baseName2="${fastq2##*/}"
baseName2="${baseName2%.*}"

final_baseName="${baseName1%_*}"
final_baseName="${final_baseName%_*}"
final_baseName="${final_baseName%_*}"


################################################################################
#Set up intermediate final names.
################################################################################

fastq1_headers=${outDir}"/${baseName1}.header_info.txt"
fastq2_headers=${outDir}"/${baseName2}.header_info.txt"
penultimate_output=${outDir}"/${final_baseName}.tmi.txt"



################################################################################
#Set up final output File.
################################################################################

final_output=${outDir}"/${final_baseName}.all_info.txt"



################################################################################
#Extract the info from the files.
#For each of the input paired FASTQ files, we extract every 4th line,
#starting with the first with the command:
#awk '(NR+3) % 4 == 0'
#We then feed that into a few commands that just replace some strings with other strings.
#The result is a tab-delimited text file which has
#the sequence name, the barcode, and the index in its 3 columns.
################################################################################

cat ${fastq1} | awk '(NR+3) % 4 == 0' | sed s/\ 1:N:0:/\\t/g | sed s/+/\\t/g > ${fastq1_headers}
cat ${fastq2} | awk '(NR+3) % 4 == 0' | sed s/\ 2:N:0:/\\t/g | sed s/+/\\t/g > ${fastq2_headers}

################################################################################
#We next paste the two header info files together.
#This produces a file with the following information:
#Header, Barcode1, Index1, Header, Barcode2, Index2
#We remove the intermediate files to avoid taking up unecessary space.
################################################################################

paste ${fastq1_headers} ${fastq2_headers} > ${penultimate_output}

rm ${fastq1_headers} ${fastq2_headers}

################################################################################
#We would like to remove that 4th column, since the header information is repeated.
#We remove that penultimate file.
################################################################################

awk '{print $1,$2,$3,$5,$6}' ${penultimate_output} > ${final_output}
rm ${penultimate_output}




#
