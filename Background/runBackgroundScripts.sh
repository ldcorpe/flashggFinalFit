#!/bin/bash


#bash variables
FILE="";
EXT="auto"; #extensiom for all folders and files created by this script
PROCS="ggh"
#CATS="UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"
#CATS="UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag"
CATS="UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,VHHadronicTag,VHTightTag,VHLooseTag"
SCALES="HighR9EE,LowR9EE,HighR9EB,LowR9EB"
SMEARS="HighR9EE,LowR9EE,HighR9EBRho,LowR9EBRho,HighR9EBPhi,LowR9EBPhi"
FTESTONLY=0
PSEUDODATAONLY=0
PSEUDODATADAT=""
SIGFILE=""
BKGPLOTSONLY=0
SEED=0

usage(){
	echo "The script runs background scripts:"
		echo "options:"

echo "-h|--help)"
echo "-i|--inputFile)"
echo "-p|--procs ) (default= ggh)"
echo "-f|--flashggCats) (defualt= UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag)"
echo "--ext)  (default= auto)"
echo "--fTestOnly) "
echo "--pseudoDataOnly) "
echo "--pseudoDataDat)"
echo "--sigFile) "
echo "--bkgPlotsOnly)"
echo "--seed) for pseudodata random number gen seed (default $SEED)"
}


#------------------------------ parsing


# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o hi:p:f: -l help,inputFile:,procs:,flashggCats:,ext:,fTestOnly,pseudoDataOnly,bkgPlotsOnly,pseudoDataDat:,sigFile:,seed: -- "$@")
then
# something went wrong, getopt will put out an error message for us
exit 1
fi
set -- $options

while [ $# -gt 0 ]
do
case $1 in
-h|--help) usage; exit 0;;
-i|--inputFile) FILE=$2; shift ;;
-p|--procs) PROCS=$2; shift ;;
-f|--flashggCats) CATS=$2; shift ;;
--ext) EXT=$2; echo "test" ; shift ;;
--fTestOnly) FTESTONLY=1; echo "ftest" ;;
--pseudoDataOnly) PSEUDODATAONLY=1;;
--pseudoDataDat) PSEUDODATADAT=$2; shift;;
--sigFile) SIGFILE=$2; shift;;
--bkgPlotsOnly) BKGPLOTSONLY=1;;
--seed) SEED=$2; shift;;

(--) shift; break;;
(-*) usage; echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
(*) break;;
esac
shift
done


OUTDIR="outdir_$EXT"
echo "[INFO] outdir is $OUTDIR" 
#if [ "$FILE" == "" ];then
#	echo "ERROR, input file (--inputFile or -i) is mandatory!"
#	exit 0
#fi

mkdir -p $OUTDIR

if [ $FTESTONLY == 0 -a $PSEUDODATAONLY == 0 -a $BKGPLOTSONLY == 0 ]; then
#IF not particular script specified, run all!
FTESTONLY=1
PSEUDODATAONLY=1
BKGPLOTSONLY=1
fi

####################################################
################## PSEUDODATAONLY ###################
####################################################

if [ $PSEUDODATAONLY == 1 ]; then

mkdir -p $OUTDIR/pseudoData

echo "--------------------------------------"
echo "Running Pseudodata"
echo "--> Create fake data by fitting simulations, throwing toys and adding datasets"
echo "--------------------------------------"

./bin/pseudodataMaker -i $PSEUDODATADAT --pseudodata 1 --plotdir $OUTDIR/pseudoData -f $CATS --seed $SEED

FILE=pseudoWS.root

fi

####################################################
################## F-TEST ###################
####################################################
if [ $FTESTONLY == 1 ]; then

echo "--------------------------------------"
echo "Running Background F-Test"
echo "-->Greate background model"
echo "--------------------------------------"

./bin/fTest -i $FILE --saveMultiPdf CMS-HGG_multipdf_$EXT.root  -D $OUTDIR/bkgfTest -f $CATS 

fi

####################################################
################### BKGPLOTS ###################
####################################################

if [ $BKGPLOTSONLY == 1 ]; then
echo "--------------------------------------"
echo "-->Create Background Validation plots"
echo "--------------------------------------"

if [ "$SIGFILE" != "" ]; then
SIG="-s $SIGFILE"
fi
echo "./scripts/subBkgPlots.py -b CMS-HGG_multipdf_$EXT.root -d $OUTDIR/bkgPlots -S 13 --isMultiPdf --useBinnedData --runLocal -c 9 $SIG"
./scripts/subBkgPlots.py -b CMS-HGG_multipdf_$EXT.root -d $OUTDIR/bkgPlots -S 13 --isMultiPdf --useBinnedData  --doBands --runLocal -c 9 --massStep 2 $SIG -L 120 -H 130 -f $CATS #for now

fi


if [ $USER == "lcorpe" ]; then
cp -r $OUTDIR ~/www/.
cp ~lcorpe/public/index.php ~/www/$OUTDIR/pseudoData/.
cp ~lcorpe/public/index.php ~/www/$OUTDIR/bkgPlots/.
cp ~lcorpe/public/index.php ~/www/$OUTDIR/bkgfTest/.

echo "plots available at: "
echo "https://lcorpe.web.cern.ch/lcorpe/$OUTDIR"

fi