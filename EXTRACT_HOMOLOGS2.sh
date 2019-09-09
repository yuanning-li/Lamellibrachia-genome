#!/bin/env bash

### PROGRAM VARIABLES ###
DIAMOND=`which diamond`
CDHIT=`which cd-hit`
SELECT_CONTIGS=`find ~ select_contigs.pl | egrep -m 1 "select_contigs.pl"`
INTERPROSCAN=`which interproscan.sh`
#########################

### MANDATORY INPUT FILES ###
NAMING_TEMPLATE=`echo "NULL"`
HOMOLOG_TARGETS=`echo "NULL"`
DATABASE=`echo "NULL"`
#############################

### OPTIONAL ARGUMENTS ###
SEQUENCES=`pwd | sed 's/$/\//'`
THREADS=`echo "1"`
KEEP_INTERMEDIATE=`echo "FALSE"`
OUTPUT_DIRECTORY=`echo "Extract_homologs_output"`
##########################

### FLAG ENTRY ###
while getopts ":d:kn:o:s:S:t:T:h" opt; do
	case $opt in
		d)
			DATABASE=$OPTARG
			;;
		k)
			KEEP_INTERMEDIATE=`echo "TRUE"`	
			;;
		h)
			printf "\n\tUSAGE\n\n"
			printf "\t-d <str> \tPath to diamond database file [MANDATORY]\n\n"
                        printf "\t-h\t\tHelp information and USAGE\n\n"
                        printf "\t-k\t\tKeep temporary files !WARNING: Disk-space intensive! [OPTIONAL]\n\n"
                        printf "\t-n <str> \tNaming template [MANDATORY]\n\n"
                        printf "\t-o <str> \tOutput directory name [OPTIONAL, DEFAULT = ./Extract_homologs_output]\n\n"
                        printf "\t-s <str> \tDirectory containing peptide sequence files [OPTIONAL, DEFAULT = pwd]\n\n"
                        printf "\t-S <str> \tPath to select_contigs.pl [OPTIONAL, DEFAULT = pwd]\n\n"
                        printf "\t-t <int> \tNumber of threads [DEFAULT = 1]\n\n"
                        printf "\t-T <str> \tList of target SwissProt protein homologs [MANDATORY]\n\n"
			exit 1
			;;
		n)	
			NAMING_TEMPLATE=$OPTARG
			;;
		o)
			OUTPUT_DIRECTORY=$OPTARG
			;;
		s)
			SEQUENCES=$OPTARG
			;;
		S)
			SELECT_CONTIGS=$OPTARG
			;;
		t)
			THREADS=$OPTARG
			;;
		T)
			HOMOLOG_TARGETS=$OPTARG
			;;
		\?) 
			echo "Invalid option: -$OPTARG"
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." 
			exit 1
			;;
	esac
done
######################

### CHECKING FOR MANDATORY FLAGS AND PROPER INPUT ###
if [[ $DATABASE == "NULL" ]]; then
	printf  "\n\tERROR: Path to diamond database [-d] not entered [MANDATORY].\n"
	ERROR=`echo "TRUE"`
fi 

if ! [[ -s "${DATABASE}" ]]; then
	printf  "\n\tERROR: Diamond database either does not exist or is empty.\n"
        ERROR=`echo "TRUE"`
fi

if [[ $NAMING_TEMPLATE == "NULL" ]]; then
	printf "\n\tERROR: Naming template [-n] not entered [MANDATORY].\n"
	ERROR=`echo "TRUE"`
fi 

if ! [[ -s "${NAMING_TEMPLATE}" ]]; then
        printf  "\n\tERROR: Naming template either does not exist or is empty.\n"
        ERROR=`echo "TRUE"`
fi

if [[ $HOMOLOG_TARGETS == "NULL" ]]; then
	printf "\n\tERROR: List of target SwissProt protein homologs [-T] not entered [MANDATORY].\n"
	ERROR=`echo "TRUE"`
fi

if ! [[ -s "${HOMOLOG_TARGETS}" ]]; then
        printf  "\n\tERROR: List of target SwissProt accessions either does not exist or is empty.\n"
        ERROR=`echo "TRUE"`
fi

if ! [[ -x "${SELECT_CONTIGS}" ]]; then
        printf "\n\tERROR: select_contigs.pl could not be found or is not executable.\n\t Add to path, specify location with [-S], or make executable with chmod.\n"
        ERROR=`echo "TRUE"`
fi

FASTA_COUNT=`ls -l $SEQUENCES | egrep -c "*.fasta"`
if [[ $FASTA_COUNT == "0" ]]; then
	printf "\n\tERROR: No .fasta files were found in working directory.\n\tUse [-s] to specify a directory containing .fasta files.\n\t"
	ERROR=`echo "TRUE"`
fi

if [[ -d "${OUTPUT_DIRECTORY}" ]]; then
	printf "\n\tERROR: Output directory exists. Please rename and try again.\n\t"
	ERROR=`echo "TRUE"`
fi

if [[ $ERROR == "TRUE" ]]; then
	printf "\n\tProvide [-h] flag for USAGE\n"
	printf "\tExiting...\n\n"
	exit 1
fi

####################################

### CREATE THE WORKING DIRECTORY ###
mkdir $OUTPUT_DIRECTORY
cd $OUTPUT_DIRECTORY
DIRECTORY_PATH=`pwd` #ABSOLUTE PATH VERSION OF $OUTPUT_DIRECTORY
cd ..
mkdir $DIRECTORY_PATH/Putative_homologs #MAKE FINAL OUTPUT DIRECTORY
####################################

### RENAMING FUNCTION ###
function rename_header { #OUTPUTS A NEW FILE CALLED ${NEW_ID}_${GENE}.fasta
    #PARAMETERS: 1 = RENAMING FILE (INPUT $1 ABOVE), 2 = HOMOLOGUE NAME (LINE IN INPUT $2 ABOVE), 3 = HOMOLOGUE FASTA

        TEMPLATE=`echo ${1}` #RENAMING TEMPLATE
        GENE=`echo ${2}` #HOMOLOGUE TO BE RENAMED AS
        F=`echo ${3}` #FASTA FILE CONTAINING HOMOLOGOUS SEQUENCES

        NEW_ID=`grep "${ID}" ${TEMPLATE} | awk '{print $2}'` # $ID COMES FROM ID VARIABLE CREATED IN MAIN

        awk '/^>/{print ">NEW_ID_GENE_" ++i; next}{print}' < ${F} > $DIRECTORY_PATH/${NEW_ID}_${GENE}.fasta
        sed -i "/>/s/NEW_ID/${NEW_ID}/" $DIRECTORY_PATH/${NEW_ID}_${GENE}.fasta
        sed -i "/>/s/GENE/${GENE}/" $DIRECTORY_PATH/${NEW_ID}_${GENE}.fasta
}
#########################

### LOOP CHECKS RENAMING TEMPLATE FOR LINES WITHOUT A PAIRED OLD-NEW NAMES ###
while read NAME
do
        COLUMN_NUMBER=`echo "${NAME}" | wc -w`
        if [ "${COLUMN_NUMBER}" -ne 2 ]; then
                printf "\n\tRenaming template contains a line without a "old-new" name pair\n\tExiting...\n\n"
                exit 1
        fi
done < $NAMING_TEMPLATE
##############################################################################

### INITIATING THE SCRIPT LOG ###
TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'` #FOR TRACKING PROGRESS OF SCRIPT
printf "\n####VARIABLE INPUT CHECK#### \n" > $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "Database set as:\t$DATABASE\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "Sequences found in:\t$SEQUENCES\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "Output directory:\t$DIRECTORY_PATH\n" >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "Diamond:\t\t$DIAMOND\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "cd-hit:\t\t\t$CDHIT\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "select_contigs.pl:\t$SELECT_CONTIGS\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "interproscan.sh:\t$INTERPROSCAN\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "Threads:\t\t$THREADS\n" >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "############################\n"  >> $DIRECTORY_PATH/SCRIPT_LOG.txt
printf "\n${TIME_i} - Beginning analysis....\n" >> $DIRECTORY_PATH/SCRIPT_LOG.txt
#################################

### INITIAL CLUSTERING ###
for FILE in $SEQUENCES*fasta #BEGINS THE LOOP THROUGH ALL THE FILES IN THE SPECIFIED DATASET DIRECTORY
do
	PEPTIDE_FILE=`echo "${FILE}" | awk -F"/" '{print $NF}'` #VARIABLE EXLCUDES THE PATH-TO-THE FILE WHICH IS INCLUDED IN THE $FILE VARIABLE
        ID=`echo "${PEPTIDE_FILE}" | awk -F "." '{print $1}'` #CREATE AN ID FOR EACH PEPTIDE FILE - CORRESPONDS TO COLUMN ONE IN INPUT 1

        TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'` #FOR TRACKING PROGRESS OF SCRIPT
        echo "${TIME_i} - Processing ${ID}. Beginning clustering." >> $DIRECTORY_PATH/SCRIPT_LOG.txt

        sed '/>/s/ .*//' $FILE > $DIRECTORY_PATH/${PEPTIDE_FILE}.clean_headers #CLEAN HEADERS
        $CDHIT -i $DIRECTORY_PATH/${PEPTIDE_FILE}.clean_headers -o $DIRECTORY_PATH/${PEPTIDE_FILE}.clean_headers.cd-hit100 -c 1.00 #CLUSTER SEQUENCES BY 100% IDENTITY TO REMOVE FRAGMENTS OF LARGER PROTEINS
	
	FILE=`ls $DIRECTORY_PATH | egrep "${PEPTIDE_FILE}.clean_headers.cd-hit100$"` #CHANGES FILE VARIABLE TO THE CLUSTERED DATASET

	if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
		rm -f $DIRECTORY_PATH/${PEPTIDE_FILE}.clean_headers.cd-hit100.clstr
                rm -f $DIRECTORY_PATH/${PEPTIDE_FILE}.clean_headers
        fi
        
	TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
        echo "${TIME_i} - Querying ${ID} against SwissProt database with Diamond BLASTp." >> $DIRECTORY_PATH/SCRIPT_LOG.txt

# FILTER .PEP FILE BY PRIMARY SEQUENCE HOMOLOGY (BLAST); MUST CHANGE PATH TO SWISSPROT BLASTDB #
	
	${DIAMOND} blastp --db ${DATABASE} --query ${DIRECTORY_PATH}/${FILE} --outfmt 6 --out ${DIRECTORY_PATH}/${ID}.blastp_vs_swissprot.outfmt6 --threads ${THREADS} #BLAST QUERY FILE AGAINST SWISSPROT DATABASE"
	
        TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
        echo "${TIME_i} - Finished blastp on ${ID}. Proceeding with Diamond Blastp output filtering" >> $OUTPUT_DIRECTORY/SCRIPT_LOG.txt

        sort -k1,1 -k11,11g $DIRECTORY_PATH/${ID}.blastp_vs_swissprot.outfmt6 | sort -u -k1,1 > $DIRECTORY_PATH/${ID}.besthit_vs_swissprot.outfmt6 #PULL BEST HITS FROM THE BLAST OUTPUT
        awk '{print $1, $2, $11}' $DIRECTORY_PATH/${ID}.besthit_vs_swissprot.outfmt6 > $DIRECTORY_PATH/${ID}.seq_hit_evalue.blast_out #REDUCE BLAST OUTPUT TO THREE COLUMNS - QUERY ID, SUBJECT ID, AND E-VALUE
        
	if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
		rm -f $DIRECTORY_PATH/${ID}.blastp_vs_swissprot.outfmt6
	        rm -f $DIRECTORY_PATH/${ID}.besthit_vs_swissprot.outfmt6
	fi

# GO THROUGH THE LIST OF GENES OF INTEREST ($2) AND PULL PUTATIVE HOMOLOGUES AND RENAME FOR EACH SEQUENCE; MUST CHANGE PATH TO SELECT_CONTIGS.PL #

        TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
        echo "${TIME_i} - Creating fasta file of homologues of interest for ${ID}" >> $DIRECTORY_PATH/SCRIPT_LOG.txt

        while read PROT #PARSE, LINE BY LINE, A LIST OF PROTEINS HOMOLOGUES
        do
                if [ "${PROT}" == "" ]; then
                        echo "Empty line skipped"
                        continue
                fi
		grep "|${PROT}" $DIRECTORY_PATH/${ID}.seq_hit_evalue.blast_out | awk '{print $1}' > $DIRECTORY_PATH/${ID}.${PROT}.tmp_before_renaming #PULL QUERY SEQ HEADERS WITH BEST HIT TO PROTEINS OF INTEREST
                
		if [[ ! -s $DIRECTORY_PATH/${ID}.${PROT}.tmp_before_renaming ]]; then #LOOP CHECKS FOR EMPTY FILES; IF FILE IS EMPTY, PRINT USEFUL INFORMATION TO LOG FILE AND ITERATE ENCLOSING LOOP
			echo "${TIME_i} - No best hits for $PROT in $ID" >> $DIRECTORY_PATH/SCRIPT_LOG.txt
			continue 1
		fi
		
		$SELECT_CONTIGS -n $DIRECTORY_PATH/${ID}.${PROT}.tmp_before_renaming $DIRECTORY_PATH/${FILE} $DIRECTORY_PATH/${ID}.${PROT}.fasta.to_be_renamed #PULL CORRESPONDING SEQUENCES FROM .PEP FILE
                
		if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
			rm -f $DIRECTORY_PATH/${ID}.${PROT}.tmp_before_renaming #CLEAN-UP
       		fi

		rename_header ${NAMING_TEMPLATE} ${PROT} $DIRECTORY_PATH/${ID}.${PROT}.fasta.to_be_renamed #RENAME SEQUENCES AND FILE TO INDICATE PROTEIN HOMOLOGY; ALSO CREATES $NEW_ID VARIABLE
                
		if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
			rm -f $DIRECTORY_PATH/${ID}.${PROT}.fasta.to_be_renamed #CLEAN-UP
		fi
        
	done < $HOMOLOG_TARGETS

	if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
	        rm -f $DIRECTORY_PATH/${ID}.seq_hit_evalue.blast_out
	fi

# CONSOLIDATE HOMOLOGUE FILES INTO ONE PER TAXON #

        TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
        echo "${TIME_i} - Renamed ${ID} to ${NEW_ID} and consolidating sequences" >> $DIRECTORY_PATH/SCRIPT_LOG.txt

        ls $DIRECTORY_PATH | grep "${NEW_ID}.*fasta$" > $DIRECTORY_PATH/${NEW_ID}_homologs.txt #ASSEMBLE A FILE CONTAINING THE FILE NAMES FOR ALL PROTEINS PULLED FROM PROTEIN DATASET IN BLOCK ABOVE
        
	while read HOMOLOG_OUTPUT_FILE
        do
                cat $DIRECTORY_PATH/$HOMOLOG_OUTPUT_FILE >> $DIRECTORY_PATH/Putative_homologs/${NEW_ID}_putative_homologs.fasta #CONSOLIDATE ALL PUTATIVE HOMOLOGS INTO ONE FILE PER TAXON
        done < $DIRECTORY_PATH/${NEW_ID}_homologs.txt

	NUMBER_OF_SEQUENCES=`grep -c '>' $DIRECTORY_PATH/Putative_homologs/${NEW_ID}_putative_homologs.fasta`
	echo "${TIME_i} - Identified $NUMBER_OF_SEQUENCES putative homologs in $NEW_ID" >> $DIRECTORY_PATH/SCRIPT_LOG.txt

	if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
	        rm -f $DIRECTORY_PATH/${NEW_ID}*.fasta #CLEAN-UP
        	rm -f $DIRECTORY_PATH/${NEW_ID}_homologs.txt #CLEAN-UP
        	rm -f $DIRECTORY_PATH/${FILE}
		find $DIRECTORY_PATH -empty -delete
	fi
done

# NOW WORKING WITH A CONSOLIDATED ALL-TAXA DATASET #

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Finished processing all datasets. Consolidating and annotating with Interproscan." >> $DIRECTORY_PATH/SCRIPT_LOG.txt

cat $DIRECTORY_PATH/Putative_homologs/*.fasta >> $DIRECTORY_PATH/Putative_homologs/Putative_homologs_from_all_taxa.fasta #CONCATONATE ALL FASTA FILES INTO ONE

# FINAL BLAST FOR CONSISE PRIMARY SEQUENCE HOMOLOGY FOR TARGET HOMOLOGS #

${DIAMOND} blastp --db ${DATABASE} --query ${DIRECTORY_PATH}/Putative_homologs/Putative_homologs_from_all_taxa.fasta --outfmt 6 \
--out ${DIRECTORY_PATH}/Putative_homologs/Putative_homologs_from_all_taxa.blastp_vs_swissprot.outfmt6 --threads ${THREADS}
sort -k1,1 -k11,11g $DIRECTORY_PATH/Putative_homologs/Putative_homologs_from_all_taxa.blastp_vs_swissprot.outfmt6 | sort -u -k1,1 > $DIRECTORY_PATH/Putative_homologs/Putative_homologs_from_all_taxa.besthit_vs_swissprot.outfmt6 #PULL BEST HITS FROM THE BLAST OUTPUT

# ANNOTATED WITH INTERPROSCAN - PFAM AND SMART #

mkdir $DIRECTORY_PATH/tmp #TMP DIRECTORY FOR INTERPRO TMP FILES
sed 's/*//g' $DIRECTORY_PATH/Putative_homologs/Putative_homologs_from_all_taxa.fasta > $DIRECTORY_PATH/Putative_homologs/Putative_homologs_from_all_taxa.interpro #INTERPRO REQUIRES NO "*", THIS STEP REMOVES THEM FROM THE FASTA FILE
$INTERPROSCAN -appl SMART -appl Pfam -cpu ${THREADS} -T $DIRECTORY_PATH/tmp -d $DIRECTORY_PATH/Putative_homologs/ -i $DIRECTORY_PATH/Putative_homologs/Putative_homologs_from_all_taxa.interpro
rmdir $DIRECTORY_PATH/tmp

# FINAL CLEANUP #
if [[ $KEEP_INTERMEDIATE == "FALSE" ]]; then
	rm -f Putative_homologs_from_all_taxa.interpro
	rm -f $DIRECTORY_PATH/$FILE
fi

TIME_i=`date | awk '{print $2, $3, $6, $4, $5}'`
echo "${TIME_i} - Completed pipeline. Good Job." >> $DIRECTORY_PATH/SCRIPT_LOG.txt

