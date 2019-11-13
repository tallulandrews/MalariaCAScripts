#!/bin/bash
# Arguments:
# $1 = Annotation GTF file
# $2 = number of threads to run on
# $3 = input BAM to run on (only required if number of threads > 0)
# $4 = workingdir
# $5 = outputdir
# This runs fast & efficiently, not on cluster took < 10 minutes to count one of the dedupped merged files.
# But should run this on the complete annotations after cufflinks de novo transcript assembly.
# To do: modify this so multi/not and single/paired are argument options

featureCOUNT=/nfs/users/nfs_t/ta6/RNASeqPipeline/software/subread-1.4.6-p2-Linux-x86_64/bin/featureCounts
ANNOTATIONgtf=$1
NUMTHREADS=$2
INPUTDIR=$3
WORKINGDIR=$4/$LSB_JOBINDEX
OUTDIR=$5
FILEStoMAP=($INPUTDIR/*.bam)
ARRAYINDEX=$(($LSB_JOBINDEX-1))
INPUTBAM=${FILEStoMAP[$ARRAYINDEX]}
echo "Inputfile: $INPUTBAM"



if [ ! -f $featureCOUNT ] ; then
  echo "Sorry featureCounts not available"
  exit 1
fi

if [ -z $ANNOTATIONgtf ] || [ ! -f $ANNOTATIONgtf ] ; then
  echo "Please provide an annotation GTF file (ARG 1/3)"
  exit 1
fi

if [ -z $NUMTHREADS ] ; then
  echo "Please set number of threads to run on (ARG 2/3)"
  exit 1
fi

if [ $NUMTHREADS -gt 0 ] ; then
  if [ -z $INPUTBAM ] || [ ! -f $INPUTBAM ] ; then
    echo "$INPUTBAM does not exist. Please provide existing sorted BAM file (ARG 3/3)"
    exit 1
  fi
else
  echo "Error: number of threads must be > 0"
  exit 1
fi

mkdir -p $WORKINGDIR
cd $WORKINGDIR

OUTPUTFILE=$(basename "${INPUTBAM%.bam}.fragmentcounts")
$featureCOUNT -O -M --minReadOverlap 70 -t exon -f -T $NUMTHREADS -a $ANNOTATIONgtf -o $OUTDIR/$OUTPUTFILE $INPUTBAM 
rm temp*
cd ..
rmdir $WORKINGDIR
