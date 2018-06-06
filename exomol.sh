#!/bin/bash

# This script downloads and unpacks  the *.states, *.pf and, *.trans
# and *.def files from wwww.exomol.com.
# And it generates the information for the ISO.h file for heliosk.
# February 2018
# Author: Simon Grimm

#run with "bash exomol.sh id" where id is the number of the molecule


m=$1				#Molecule id
PrintISO=1			#when set to 1, then print the code fot the ISO.h file
DownloadFiles=1

echo "molecule "$m 

if [ $m -eq 1 ]
then
  #1 H2O
  M="1H2-16O__BT2"
  P="H2O/1H2-16O/BT2"
  s=-1				#file range
  ntcol=3			#columns in transition files
  npfcol=3			#columns in partition function file
fi
if [ $m -eq 101 ]
then
  #the number of transition files is from, correct later
  #1 H2O
  M="1H2-16O__POKAZATEL"
  P="H2O/1H2-16O/POKAZATEL/"
  s=100				#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 5 ]
then
  #5 CO
  M="12C-16O__Li2015"
  P="CO/12C-16O/Li2015"
  s=22000			#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi


if [ $m -eq 6 ]
then
  #6 CH4
  M="12C-1H4__YT10to10"
  P="CH4/12C-1H4/YT10to10"
  s=100				#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 8 ]
then
  #8 NO
  M="14N-16O__NOname"
  P="NO/14N-16O/NOname"
  s=40000				#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 9 ]
then
  #9 SO2
  M="32S-16O2__ExoAmes"
  P="SO2/32S-16O2/ExoAmes"
  s=100				#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi


if [ $m -eq 11 ]
then
  #11 NH3
  M="14N-1H3__BYTe"
  P="NH3/14N-1H3/BYTe"
  s=100				#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 23 ]
then
  #23 HCN
  M="1H-12C-14N__Harris"
  P="HCN/1H-12C-14N/Harris"
  s=17586			#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 28 ]
then
  #28 PH3
  M="31P-1H3__SAlTY"
  P="PH3/31P-1H3/SAlTY"
  s=100				#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 31 ]
then
  #31 H2S
  M="1H2-32S__AYT2"
  P="H2S/1H2-32S/AYT2"
  s=1000			#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 80 ]
then
  #80 VO
  M="51V-16O__VOMYT"
  P="VO/51V-16O/VOMYT"
  s=5000			#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

#81 TiO

if [ $m -eq 82 ]
then
  #82 FeH
  M="56Fe-1H__Yueqi"
  P="FeH/56Fe-1H/Yueqi"
  s=7476			#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 83 ]
then
  #83 AlO
  M="27Al-16O__ATP"
  P="AlO/27Al-16O/ATP"
  s=35000			#file range
  ntcol=3			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 84 ]
then
  #84 SiO
  M="28Si-16O__EBJT"
  P="SiO/28Si-16O/EBJT"
  s=6050			#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 85 ]
then
  #85 CaO
  M="40Ca-16O__VBATHY"
  P="CaO/40Ca-16O/VBATHY"
  s=25000			#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 86 ]
then
  #86 SiH
  M="28Si-1H__SiGHTLY"
  P="SiH/28Si-1H/SiGHTLY"
  s=32000			#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 87 ]
then
  #87 caH
  M="40Ca-1H__Yadin"
  P="CaH/40Ca-1H/Yadin"
  s=15278			#file range
  ntcol=4			#columns in transition files
  npfcol=2			#columns in partition function file
fi

if [ $m -eq 88 ]
then
  #88 H3+
  M="1H3_p__MiZATeP"
  P="H3_p/1H3_p/MiZATeP"
  s=25000			#file range
  ntcol=3			#columns in transition files
  npfcol=3			#columns in partition function file
fi

if [ $m -eq 89 ]
then
  #89 PO
  M="31P-16O__POPS"
  P="PO/31P-16O/POPS"
  s=12000			#file range
  ntcol=4			#columns in transition files
  npfcol=3			#columns in partition function file
fi
if [ $m -eq 90 ]
then
  #90 MgH
  M="24Mg-1H__Yadin"
  P="MgH/24Mg-1H/Yadin"
  s=11000			#file range
  ntcol=4			#columns in transition files
  npfcol=3			#columns in partition function file
fi


echo $M

if [ $DownloadFiles -eq 1 ]
then
wget http://exomol.com/db/$P/$M.states.bz2
bzip2 -d $M.states.bz2
wget http://exomol.com/db/$P/$M.pf
fi
wget http://exomol.com/db/$P/$M.def

n=`grep "No. of transition files" $M.def | cut -c-12`
mass=`grep "Isotopologue mass" $M.def | cut -c-12`
dL=`grep "Default value of Lorentzian half-width for all lines" $M.def | cut -c-12`
dn=`grep "Default value of temperature exponent for all lines" $M.def | cut -c-12`
version=`grep "Version number with format" $M.def | cut -c-12`

if [ $m -eq 101 ]
then
n=412

fi

echo $mass
echo $dL
echo $dn

for (( nu=0; nu<$n; nu++ ))
do
  echo "------Download file "$nu" from "$n" -----"

  if [ $m -eq 1 ]
  then
    jarray[0]=00000
    jarray[1]=00250
    jarray[2]=00500
    jarray[3]=00750
    jarray[4]=01000
    jarray[5]=01500
    jarray[6]=02000
    jarray[7]=02250
    jarray[8]=02750
    jarray[9]=03500
    jarray[10]=04500
    jarray[11]=05500
    jarray[12]=07000
    jarray[13]=09000
    jarray[14]=14000
    jarray[15]=20000
    jarray[16]=30000

    if [ $DownloadFiles -eq 1 ]
    then
      wget http://www.exomol.com/db/$P/$M\_\_${jarray[$nu]}-${jarray[$nu + 1]}.trans.bz2
      bzip2 -d $M\_\_${jarray[$nu]}-${jarray[$nu + 1]}.trans.bz2
    fi
    l[$nu]=`wc -l < $M\_\_${jarray[$nu]}-${jarray[$nu + 1]}.trans | awk '{print $1}'`

  else
    printf -v j "%05d" $((nu*$s))
    printf -v jj "%05d" $(($nu*$s+$s))
    if [ $n -gt 1 ]
    then
      if [ $DownloadFiles -eq 1 ]
      then
        wget http://www.exomol.com/db/$P/$M\_\_$j-$jj.trans.bz2
        bzip2 -d $M\_\_$j-$jj.trans.bz2
      fi
      l[$nu]=`wc -l < $M\_\_$j-$jj.trans | awk '{print $1}'`
    else
      if [ $DownloadFiles -eq 1 ]
      then
        wget http://www.exomol.com/db/$P/$M.trans.bz2
        bzip2 -d $M.trans.bz2
      fi
      l[$nu]=`wc -l < $M.trans | awk '{print $1}'`
    fi
  fi

  echo $nu ${l[$nu]}
done

echo download finished

if [ $PrintISO -eq 1 ]
then
  ll=`wc -l < $M.states | awk '{print $1}'`
  echo char name"[]" = \"$M\"";"
  echo -e "sprintf(m.mName, \"%s\", \"$M\");"
  echo m.defaultL = $dL";"
  echo m.defaultn = $dn";"
  echo m.nStates = $ll";"
  echo m.nFiles = $n";"
  echo m.ntcol = $ntcol";"
  echo m.npfcol = $npfcol";"
  for (( nu=0; nu<$n; nu++ ))
  do
    echo m.NL[$nu] = ${l[$nu]}";" 
  done
  echo -e 'm.NLmax = 0;'
  if [ $m -eq 1 ]
  then
    for (( nu=0; nu<$n+1; nu++ ))
    do
      echo m.fileLimit[$nu] = ${jarray[$nu]}";" 
    done

  else
    echo -e 'for(int i = 0; i < m.nFiles + 1; ++i){'
    echo  "	m.fileLimit[i] = i * "$s";" 
    echo -e '	m.NLmax = max(m.NLmax, m.NL[i]);'
    echo -e '}'
  fi

  echo -e 'sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");'

  if [ $n -gt 1 ]
  then
    echo -e 'for(int i = 0; i < m.nFiles; ++i){'
    echo -e '	sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);'
    echo -e '}'
  else
    echo -e '	sprintf(m.dataFilename[0], "%s%s.", param.path, name);'

  fi

  echo -e 'm.nISO = 1;'
  echo -e 'm.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));'
  echo -e "m.ISO[0] = (Isotopologue){XX1,  XX,  1.0,    0.0,    0,     "$mass"};"

  echo -e "version = "$version
fi

