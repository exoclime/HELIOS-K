#!/bin/bash

# This script downloads and unpacks  the *.states, *.pf and, *.trans
# and *.def files from wwww.exomol.com.
# And it generates the information for the ISO.h file for heliosk.
# November 2017
# Author: Simon Grimm


M=$1

echo "molecule "$M 

if [ $M == 6 ]
then
  #6 CH4
  n=121				#number of files
  M="12C-1H4__YT10to10"
  P="CH4/12C-1H4/YT10to10"
  s=100				#file range
  ntcol = 3
fi


if [ $M == 11 ]
then
  #11 NH3
  n=120				#number of files
  M="14N-1H3__BYTe"
  P="NH3/14N-1H3/BYTe"
  s=100				#file range
  ntcol = 3
fi

if [ $M == 23 ]
then
  #23 HCN
  n=1				#number of files
  M="1H-12C-14N__Harris"
  P="HCN/1H-12C-14N/Harris"
  s=17586			#file range
  ntcol = 4
fi

if [ $M == 31 ]
then
  #31 H2S
  n=35				#number of files
  M="1H2-32S__AYT2"
  P="H2S/1H2-32S/AYT2"
  s=1000			#file range
  ntcol = 3
fi

if [ $M == 80 ]
then
  #80 VO
  n=7				#number of files
  M="51V-16O__VOMYT"
  P="VO/51V-16O/VOMYT"
  s=5000			#file range
  ntcol = 4
fi


echo $M


wget http://exomol.com/db/$P/$M.states.bz2
bzip2 -d $M.states.bz2
wget http://exomol.com/db/$P/$M.pf
wget http://exomol.com/db/$P/$M.def


mass=`grep "Isotopologue mass" $M.def | cut -c-12`
dL=`grep "Default value of Lorentzian half-width for all lines" $M.def | cut -c-12`
dn=`grep "Default value of temperature exponent for all lines" $M.def | cut -c-12`
version=`grep "Version number with format" $M.def | cut -c-12`

echo $mass
echo $dL
echo $dn

for (( nu=0; nu<$n; nu++ ))
do
  echo "------Download file "$nu" from "$n" -----"

  printf -v j "%05d" $((nu*$s))
  printf -v jj "%05d" $(($nu*$s+$s))
  if [ $n -gt 1 ]
  then
    wget http://www.exomol.com/db/$P/$M\_\_$j-$jj.trans.bz2
    bzip2 -d $M\_\_$j-$jj.trans.bz2
  else
    wget http://www.exomol.com/db/$P/$M.trans.bz2
    bzip2 -d $M.trans.bz2
  fi
  if [ $n -gt 1 ]
  then
    l[$nu]=`wc -l < $M\_\_$j-$jj.trans | awk '{print $1}'`
  else
    l[$nu]=`wc -l < $M.trans | awk '{print $1}'`
  fi
  echo $nu ${l[$nu]}
done

if [ $true ]
then
  ll=`wc -l < $M.states | awk '{print $1}'`
  echo char name"[]" = \"$M\"";"
  echo m.defaultL = $dL";"
  echo m.defaultn = $dn";"
  echo m.nStates = $ll";"
  echo m.nFiles = $n";"
  echo m.ntcol = $ntcol";"
  for (( nu=0; nu<$n; nu++ ))
  do
    echo m.NL[$nu] = ${l[$nu]}";" 
  done
  echo -e 'm.NLmax = 0;'
  echo -e 'for(int i = 0; i < m.nFiles + 1; ++i){'
  echo	"	m.fileLimit[i] = i * "$s";"
  echo -e '	m.NLmax = max(m.NLmax, m.NL[i]);'
  echo -e '}'

  echo -e 'sprintf(qFilename, "%s%s%s", param.path, name, ".pf");'

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

