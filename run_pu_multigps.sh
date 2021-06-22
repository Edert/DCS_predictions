
NAME=$1
SET=$2
TOOL=$3
WORKER=5

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL


files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13"
  
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)

  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`
  
  echo "#DataFile	Signal/Control	FileFormat	ConditionName	ReplicateName" > design.txt
  echo "$4	signal	BAM	S1	1" >> design.txt
  echo "$5	signal	BAM	S1	2" >> design.txt
  echo "$6	control	BAM	S1	" >> design.txt
  echo "$7	signal	BAM	S2	1" >> design.txt 
  echo "$8	signal	BAM	S2	2" >> design.txt
  echo "$9	control	BAM	S2	" >> design.txt
  
  PREPDONE=`date +%s.%N`
  
  #default run #increased from 20G to 40G
  #java -Xmx40G -jar /apps/multigps_v0.74/multigps_v0.74.jar --out $OUT_NAME --threads 1 --geninfo /ssd/references/THOR/mm10/chr19.chrom.size --design design.txt --nomotifs --minfold 0.7 --diffp 1.0 >> $LOG 2>&1
  /usr/bin/time -o mem.txt -f "%K %M" java -Xms20G -Xmx60G -d64 -jar /apps/multigps_v0.74/multigps_v0.74.jar --out $OUT_NAME --threads $WORKER --geninfo /ssd/references/THOR/mm10/chr19.chrom.size --design design.txt --nomotifs --minfold 0.7 --diffp 1.0 >> $LOG 2>&1
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e $OUT_NAME/$OUT_NAME.overdisp0.15.S1vsS2.edgeR_GLM_DE.txt ]; then
    #move results
    cat $OUT_NAME/$OUT_NAME.overdisp0.15.S1vsS2.edgeR_GLM_DE.txt | sed "s/\"//g" | grep -v "^logFC" |  awk '{split($1,a,":");print a[1]"\t"a[2]-50"\t"a[2]+50"\t"$6"\t"$2 }' | sort -k1,1 -k2,2n > $OUT_NAME".bed"
  else
      #create empty file
      touch $OUT_NAME".bed"
  fi
  
  #clean up
  rm -f design.txt mem.txt
  rm -rf $OUT_NAME
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
