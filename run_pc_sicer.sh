
NAME=$1
SET=$2
TOOL=$3

mkdir -p results_peaks/$TOOL results_peaks/$TOOL/$SET
mkdir -p log/pc log/pc/$TOOL

#get effective genome size for mm10 chr19
#epic-effective --read-length=50 --nb-cpu=4 ../01_simulating/data/mm10_chr19.fasta
#Number unique 50-mers:  55405574
#Effective genome size:  0.9019072377220532

files=$(ls results_peaks/$TOOL/$SET/*/*.bed 2> /dev/null | wc -l)
if [ "$files" = "0" ]; then

  echo "running peak calling tool: $TOOL with: $1 $2 $3" #$4 $5 $6 $7 $8 $9"
  
  cd results_peaks/$TOOL/$SET/
  
  LOG="../../../log/pc/$TOOL/$SET.log"
  RESULTSDIR="../../../results_peaks/$TOOL/$SET"

  bamToBed -i $4 > S11.bed
  bamToBed -i $5 > S12.bed
  bamToBed -i $6 > IN1.bed
  bamToBed -i $7 > S21.bed
  bamToBed -i $8 > S22.bed
  bamToBed -i $9 > IN2.bed

  #sharp
  WINDOW=50
  GAP=100
  PCRUN="sharp"
  mkdir -p $PCRUN
  
  sicer -t S11.bed -c IN1.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP -cpu 1 >> $LOG 2>&1
  cat S11-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s11_peaks.bed #chr, start, stop, 1-fdr to get some weight
  #clean up
  rm -f S11-W$WINDOW-G$GAP* S11-W$WINDOW-normalized.wig 
  
  sicer -t S12.bed -c IN1.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S12-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s12_peaks.bed 
  #clean up
  rm -f S12-W$WINDOW-G$GAP* S12-W$WINDOW-normalized.wig 
  
  sicer -t S21.bed -c IN2.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S21-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s21_peaks.bed 
  #clean up
  rm -f S21-W$WINDOW-G$GAP* S21-W$WINDOW-normalized.wig 
  
  sicer -t S22.bed -c IN2.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S22-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s22_peaks.bed 
  #clean up
  rm -f S22-W$WINDOW-G$GAP* S22-W$WINDOW-normalized.wig 

  #broad
  WINDOW=100
  GAP=200
  PCRUN="broad"
  mkdir -p $PCRUN
  
  sicer -t S11.bed -c IN1.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S11-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s11_peaks.bed 
  #clean up
  rm -f S11-W$WINDOW-G$GAP* S11-W$WINDOW-normalized.wig 
  
  sicer -t S12.bed -c IN1.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S12-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s12_peaks.bed 
  #clean up
  rm -f S12-W$WINDOW-G$GAP* S12-W$WINDOW-normalized.wig 
  
  sicer -t S21.bed -c IN2.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S21-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s21_peaks.bed 
  #clean up
  rm -f S21-W$WINDOW-G$GAP* S21-W$WINDOW-normalized.wig 
  
  sicer -t S22.bed -c IN2.bed -s mm10 -rt 1 -w $WINDOW -f 200 -egf 0.9019072377220532 -fdr 0.01 -g $GAP  -cpu 1 >> $LOG 2>&1
  cat S22-W$WINDOW-G$GAP-islands-summary | awk '{ if ($8 <= 0.01) {printf "%s\t%s\t%s\t%.30f\n", $1, $2, $3, 1-$8 }}' | sort -k1,1 -k2,2n > $PCRUN/s22_peaks.bed 
  #clean up
  rm -f S22-W$WINDOW-G$GAP* S22-W$WINDOW-normalized.wig 

  #clean up
  rm -rf S11.bed S12.bed IN1.bed S21.bed S22.bed IN2.bed chr.list

else
  echo "results_peaks/$TOOL/$SET/bed already exists exiting..."
fi
