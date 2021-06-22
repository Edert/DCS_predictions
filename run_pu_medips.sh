NAME=$1
SET=$2
TOOL=$3

mkdir -p results/$TOOL results/$TOOL/$SET
mkdir -p log/$TOOL

files=$(ls results/$TOOL/$SET/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9 $10 $11 $12 $13"
  
  OUT_NAME=$(basename $NAME _sample1-rep1_mm)
  SCRIPT="script.R"
  LOG="../../../log/$TOOL/$SET.log"
  
  cd results/$TOOL/$SET/
  
  STARTTIME=`date +%s.%N`

  # in R  tmm norm
  echo "library(MEDIPS)" > $SCRIPT
  echo "library(BSgenome.Mmusculus.UCSC.mm10)" >> $SCRIPT
  echo "BSgenome='BSgenome.Mmusculus.UCSC.mm10'" >> $SCRIPT
  echo "uniq=0" >> $SCRIPT
  echo "extend=300" >> $SCRIPT
  echo "shift=0" >> $SCRIPT
  echo "ws=100" >> $SCRIPT
  echo "chr.select='chr19'" >> $SCRIPT
  echo "S11 = MEDIPS.createSet(file = '$4', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S12 = MEDIPS.createSet(file = '$5', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S1 = c(S11,S12)" >> $SCRIPT
  echo "S21 = MEDIPS.createSet(file = '$7', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S22 = MEDIPS.createSet(file = '$8' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2 = c(S21,S22)" >> $SCRIPT
  echo "S1_Input = MEDIPS.createSet(file = '$6', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2_Input = MEDIPS.createSet(file = '$9', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "CS = MEDIPS.couplingVector(pattern = 'CG', refObj = S1[[1]])" >> $SCRIPT
  echo "resultTable = MEDIPS.meth(MSet1 = S1, MSet2 = S2, CSet = CS, ISet1 = S1_Input, ISet2 = S2_Input, chr = 'chr19', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')" >> $SCRIPT
  echo "sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='results.txt', sig, quote = FALSE, sep = \"\t\",row.names = FALSE, col.names = FALSE)" >> $SCRIPT
  
  PREPDONE=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFFPREP=`echo "$PREPDONE - $STARTTIME" | bc | awk -F"." '{print}'`
  TIMEDIFF=`echo "$ENDTIME - $PREPDONE" | bc | awk -F"." '{print}'`
  echo "prep $TIMEDIFFPREP" > time.txt
  echo "1 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "1 $MEMUSAGE" >> memory.txt
  
  if [ -e results.txt ]; then
     #reformat for eval
     cat results.txt | awk '{print $1"\t"$2"\t"$3"\t"$23"\t"$20}' | sort -k1,1 -k2,2n > $OUT_NAME"_1.bed"
  else
     #create empty file
     touch $OUT_NAME"_1.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.txt mem.txt
  
  
  # in R  quantile norm
  echo "library(MEDIPS)" > $SCRIPT
  echo "library(BSgenome.Mmusculus.UCSC.mm10)" >> $SCRIPT
  echo "BSgenome='BSgenome.Mmusculus.UCSC.mm10'" >> $SCRIPT
  echo "uniq=0" >> $SCRIPT
  echo "extend=300" >> $SCRIPT
  echo "shift=0" >> $SCRIPT
  echo "ws=100" >> $SCRIPT
  echo "chr.select='chr19'" >> $SCRIPT
  echo "S11 = MEDIPS.createSet(file = '$4', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S12 = MEDIPS.createSet(file = '$5', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S1 = c(S11,S12)" >> $SCRIPT
  echo "S21 = MEDIPS.createSet(file = '$7', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S22 = MEDIPS.createSet(file = '$8' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2 = c(S21,S22)" >> $SCRIPT
  echo "S1_Input = MEDIPS.createSet(file = '$6', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2_Input = MEDIPS.createSet(file = '$9', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "CS = MEDIPS.couplingVector(pattern = 'CG', refObj = S1[[1]])" >> $SCRIPT
  echo "resultTable = MEDIPS.meth(MSet1 = S1, MSet2 = S2, CSet = CS, ISet1 = S1_Input, ISet2 = S2_Input, chr = 'chr19', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='quantile')" >> $SCRIPT
  echo "sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='results.txt', sig, quote = FALSE, sep = \"\t\",row.names = FALSE, col.names = FALSE)" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`

  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "2 $TIMEDIFF" >> time.txt 
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "2 $MEMUSAGE" >> memory.txt
  
  if [ -e results.txt ]; then
     #reformat for eval
     cat results.txt | awk '{print $1"\t"$2"\t"$3"\t"$23"\t"$20}' | sort -k1,1 -k2,2n > $OUT_NAME"_2.bed"
  else
     #create empty file
     touch $OUT_NAME"_2.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.txt mem.txt
  
  
  # in R  no norm
  echo "library(MEDIPS)" > $SCRIPT
  echo "library(BSgenome.Mmusculus.UCSC.mm10)" >> $SCRIPT
  echo "BSgenome='BSgenome.Mmusculus.UCSC.mm10'" >> $SCRIPT
  echo "uniq=0" >> $SCRIPT
  echo "extend=300" >> $SCRIPT
  echo "shift=0" >> $SCRIPT
  echo "ws=100" >> $SCRIPT
  echo "chr.select='chr19'" >> $SCRIPT
  echo "S11 = MEDIPS.createSet(file = '$4', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S12 = MEDIPS.createSet(file = '$5', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S1 = c(S11,S12)" >> $SCRIPT
  echo "S21 = MEDIPS.createSet(file = '$7', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S22 = MEDIPS.createSet(file = '$8' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2 = c(S21,S22)" >> $SCRIPT
  echo "S1_Input = MEDIPS.createSet(file = '$6', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2_Input = MEDIPS.createSet(file = '$9', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "CS = MEDIPS.couplingVector(pattern = 'CG', refObj = S1[[1]])" >> $SCRIPT
  echo "resultTable = MEDIPS.meth(MSet1 = S1, MSet2 = S2, CSet = CS, ISet1 = S1_Input, ISet2 = S2_Input, chr = 'chr19', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='none')" >> $SCRIPT
  echo "sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='results.txt', sig, quote = FALSE, sep = \"\t\",row.names = FALSE, col.names = FALSE)" >> $SCRIPT
  
  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "3 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "3 $MEMUSAGE" >> memory.txt
  
  if [ -e results.txt ]; then
     #reformat for eval
     cat results.txt | awk '{print $1"\t"$2"\t"$3"\t"$23"\t"$20}' | sort -k1,1 -k2,2n > $OUT_NAME"_3.bed"
  else
     #create empty file
     touch $OUT_NAME"_3.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.txt mem.txt
  
  
  #same with window 1000
  
  # in R  tmm norm
  echo "library(MEDIPS)" > $SCRIPT
  echo "library(BSgenome.Mmusculus.UCSC.mm10)" >> $SCRIPT
  echo "BSgenome='BSgenome.Mmusculus.UCSC.mm10'" >> $SCRIPT
  echo "uniq=0" >> $SCRIPT
  echo "extend=300" >> $SCRIPT
  echo "shift=0" >> $SCRIPT
  echo "ws=1000" >> $SCRIPT
  echo "chr.select='chr19'" >> $SCRIPT
  echo "S11 = MEDIPS.createSet(file = '$4', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S12 = MEDIPS.createSet(file = '$5', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S1 = c(S11,S12)" >> $SCRIPT
  echo "S21 = MEDIPS.createSet(file = '$7', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S22 = MEDIPS.createSet(file = '$8' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2 = c(S21,S22)" >> $SCRIPT
  echo "S1_Input = MEDIPS.createSet(file = '$6', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2_Input = MEDIPS.createSet(file = '$9', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "CS = MEDIPS.couplingVector(pattern = 'CG', refObj = S1[[1]])" >> $SCRIPT
  echo "resultTable = MEDIPS.meth(MSet1 = S1, MSet2 = S2, CSet = CS, ISet1 = S1_Input, ISet2 = S2_Input, chr = 'chr19', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='tmm')" >> $SCRIPT
  echo "sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='results.txt', sig, quote = FALSE, sep = \"\t\",row.names = FALSE, col.names = FALSE)" >> $SCRIPT

  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "4 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "4 $MEMUSAGE" >> memory.txt
  
  if [ -e results.txt ]; then
     #reformat for eval
     cat results.txt | awk '{print $1"\t"$2"\t"$3"\t"$23"\t"$20}' | sort -k1,1 -k2,2n > $OUT_NAME"_4.bed"
  else
     #create empty file
     touch $OUT_NAME"_4.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.txt mem.txt
  
  
  # in R  quantile norm
  echo "library(MEDIPS)" > $SCRIPT
  echo "library(BSgenome.Mmusculus.UCSC.mm10)" >> $SCRIPT
  echo "BSgenome='BSgenome.Mmusculus.UCSC.mm10'" >> $SCRIPT
  echo "uniq=0" >> $SCRIPT
  echo "extend=300" >> $SCRIPT
  echo "shift=0" >> $SCRIPT
  echo "ws=1000" >> $SCRIPT
  echo "chr.select='chr19'" >> $SCRIPT
  echo "S11 = MEDIPS.createSet(file = '$4', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S12 = MEDIPS.createSet(file = '$5', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S1 = c(S11,S12)" >> $SCRIPT
  echo "S21 = MEDIPS.createSet(file = '$7', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S22 = MEDIPS.createSet(file = '$8' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2 = c(S21,S22)" >> $SCRIPT
  echo "S1_Input = MEDIPS.createSet(file = '$6', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2_Input = MEDIPS.createSet(file = '$9', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "CS = MEDIPS.couplingVector(pattern = 'CG', refObj = S1[[1]])" >> $SCRIPT
  echo "resultTable = MEDIPS.meth(MSet1 = S1, MSet2 = S2, CSet = CS, ISet1 = S1_Input, ISet2 = S2_Input, chr = 'chr19', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='quantile')" >> $SCRIPT
  echo "sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='results.txt', sig, quote = FALSE, sep = \"\t\",row.names = FALSE, col.names = FALSE)" >> $SCRIPT

  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "5 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "5 $MEMUSAGE" >> memory.txt
  
  if [ -e results.txt ]; then
     #reformat for eval
     cat results.txt | awk '{print $1"\t"$2"\t"$3"\t"$23"\t"$20}' | sort -k1,1 -k2,2n > $OUT_NAME"_5.bed"
  else
     #create empty file
     touch $OUT_NAME"_5.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.txt mem.txt
  
  
  # in R  no norm
  echo "library(MEDIPS)" > $SCRIPT
  echo "library(BSgenome.Mmusculus.UCSC.mm10)" >> $SCRIPT
  echo "BSgenome='BSgenome.Mmusculus.UCSC.mm10'" >> $SCRIPT
  echo "uniq=0" >> $SCRIPT
  echo "extend=100" >> $SCRIPT
  echo "shift=0" >> $SCRIPT
  echo "ws=1000" >> $SCRIPT
  echo "chr.select='chr19'" >> $SCRIPT
  echo "S11 = MEDIPS.createSet(file = '$4', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S12 = MEDIPS.createSet(file = '$5', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S1 = c(S11,S12)" >> $SCRIPT
  echo "S21 = MEDIPS.createSet(file = '$7', BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S22 = MEDIPS.createSet(file = '$8' ,BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2 = c(S21,S22)" >> $SCRIPT
  echo "S1_Input = MEDIPS.createSet(file = '$6', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "S2_Input = MEDIPS.createSet(file = '$9', BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws, chr.select = chr.select)" >> $SCRIPT
  echo "CS = MEDIPS.couplingVector(pattern = 'CG', refObj = S1[[1]])" >> $SCRIPT
  echo "resultTable = MEDIPS.meth(MSet1 = S1, MSet2 = S2, CSet = CS, ISet1 = S1_Input, ISet2 = S2_Input, chr = 'chr19', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='none')" >> $SCRIPT
  echo "sig = MEDIPS.selectSig(results=resultTable, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)" >> $SCRIPT
  echo "options(scipen = 999)" >> $SCRIPT
  echo "write.table(file='results.txt', sig, quote = FALSE, sep = \"\t\",row.names = FALSE, col.names = FALSE)" >> $SCRIPT

  STARTTIME=`date +%s.%N`
  
  #run it...
  /usr/bin/time -o mem.txt -f "%K %M" R CMD BATCH --vanilla $SCRIPT
  
  ENDTIME=`date +%s.%N`
  TIMEDIFF=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print}'`
  echo "6 $TIMEDIFF" >> time.txt
  
  MEMUSAGE=$(sed '/non-zero status/d' mem.txt )
  echo "6 $MEMUSAGE" >> memory.txt
  
  if [ -e results.txt ]; then
     #reformat for eval
     cat results.txt | awk '{print $1"\t"$2"\t"$3"\t"$23"\t"$20}' | sort -k1,1 -k2,2n > $OUT_NAME"_6.bed"
  else
     #create empty file
     touch $OUT_NAME"_6.bed"
  fi
  
  #save log     
  cat script.Rout >> $LOG
  #clean up
  rm -f $SCRIPT script.Rout results.txt mem.txt
  
else
  echo "results/$TOOL/$SET/bed already exists exiting..."
fi
