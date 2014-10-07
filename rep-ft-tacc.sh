#!/bin/bash

BUFFSIZE=100
OUTN=$1/estimatedgenetrees


RES=()

cp $OUTN $OUTN.back

DONE=$(cat $OUTN|wc -l)
ALL=$(ls $1/*.fas|wc -l)

exec 5>>$OUTN

for i in `seq $((DONE + 1)) $ALL`; do
  x=$1/`printf '%04d' $i`.fas
  echo $i

  r=$($HOME/bin/fasttree -nt -gtr -quiet -nopr -gamma $x)
  RES+=($r)

  if [[ ${#RES[@]} == $BUFFSIZE ]]; then
    ( IFS=$'\n'; echo "${RES[*]}" )>&5
    RES=()
  fi
done
 
( IFS=$'\n'; echo "${RES[*]}" )>&5

exec 5<&-
