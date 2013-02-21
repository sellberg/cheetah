#!/usr/bin/awk -f

BEGIN {
  default=-1
  header=1    #print header
  run="NORUN?"#run number
  nxtc=-1     #Number of XTC files
  sizextc=-1  #Total size of XTC files
  darksub=-1  #Dark Cal Subtraction
  nfr=-1      #Number of Frames Processed
  nhitf=-1    #number of hits
  hitf=-1     #hitfinder
  hitfADC=-1  #hitfinderADC
  hitfNAT=-1  #hitfinderNAT
  nicef=-1    #icefinder
  icefADC=-1  #icefinderADC
  icefNAT=-1  #icefinderNAT
  hitfr=-1    #hitrate
  icefr=-1    #hitrate ice

  if (length(comment) == 0){comment "no comment"}
} 
{
  if ($0~/r0/)            {                         run=$1     }
  if ($0~/Number of XTC/) {                         nxtc=$5    }
  if ($0~/XTC total size/){                         sizextc=$5 }
  if ($0~/useDarkcalSub/) {gsub(/useDarkcal.*=/,"");darksub=$0 }
  if ($0~/Frames pro/)    {                         nfr=$3     }
  if ($0~/hitfinder=/)    {gsub(/hitfinder=/,"")   ;hitf=$0    }
  if ($0~/hitfinderADC=/) {gsub(/hitfinderADC=/,"");hitfADC=$0 }
  if ($0~/hitfinderNAT=/) {gsub(/hitfinderNAT=/,"");hitfNAT=$0 }
  if ($0~/Numb.* hits: /) {                         nhitf=$4   }
  if ($0~/Aver.*hit r/)   {                         hitfr=$4   }
  if ($0~/nFr.*ice powd/) {                         nicef=$6   }
  if ($0~/icefinder=/)    {gsub(/icefinder=/,"")   ;icef=$0    }
  if ($0~/icefinderADC=/) {gsub(/icefinderADC=/,"");icefADC=$0 }
  if ($0~/icefinderNAT=/) {gsub(/icefinderNAT=/,"");icefNAT=$0 }
}
END {
  if (nfr > 0 ) {  icefr=(nicef/nfr)*100 }
  format1="%8s ,"
  format2="%8.2f ,"
  format3="%8d ,"
  format4="%s "
  format4b="\"%s\""

  runnum = substr(run,2,4)
#  print runnum

  printf(format1,"runnum");
  printf(format1,"nxtc");
  printf(format1,"sizextc");
  printf(format1,"nfr");
  printf(format1,"darksub")
  printf(format1,"hitf")
  printf(format1,"nhitf")
  printf(format1,"hitfADC")
  printf(format1,"hitfNAT")
  printf(format1,"hitfr")
  printf(format1,"icef ")
  printf(format1,"nicef ")
  printf(format1,"icefADC")
  printf(format1,"icefNAT")
  printf(format1,"icefr")
  printf(format4,"comment")

  printf("\n")

#  printf(format1,substr(FILENAME,0,5));
  printf(format3,runnum);
  printf(format3,nxtc);
  printf(format3,sizextc)
  printf(format3,nfr);
  printf(format3,darksub)
  printf(format3,hitf)
  printf(format3,nhitf)
  printf(format3,hitfADC)
  printf(format3,hitfNAT)
  printf(format2,hitfr)
  printf(format3,icef )
  printf(format3,nicef )
  printf(format3,icefADC)
  printf(format3,icefNAT)
  printf(format2,icefr)
  printf(format4b,comment)

  printf("\n")
}

