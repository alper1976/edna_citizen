
######################################
### Fish eDNA library 221202_M01334.Project_Kvalheim-Amplicon1-2022-11-18

### cutadapt
module load cutadapt/4.2-GCCcore-11.3.0

RAWDIR=/cluster/projects/nn9745k/01_raw_data
RESDIR=/cluster/projects/nn9745k/02_results

# transfer files to scratch

mkdir $RESDIR/26_lone
mkdir $RESDIR/26_lone/AdaptersRemoved
mkdir $RESDIR/26_lone/figs


##############################################
### run cutadapt to remove primer sequences
##############################################

cd $RAWDIR/26_lone/221202_M01334.Project_Kvalheim-Amplicon1-2022-11-18


# RUN fish

for f in *.fastq.gz
    do
        FILE=${f#$DIRS}
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            cutadapt -g GTCGGTAAAACTCGTGCCAGC -G CATAGTGGGGTATCTAATCCCAGTTTG -o $RESDIR/26_lone/AdaptersRemoved/$FILE -p $RESDIR/26_lone/AdaptersRemoved/${FILE//R1/R2} --discard-untrimmed --minimum-length 35 $FILE ${FILE//R1/R2}

        fi
    done

cd $RESDIR/26_lone/AdaptersRemoved



