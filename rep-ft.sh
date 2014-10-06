#!/bin/bash

RES=()
OUTN=()

DONE=$(ls $1/*.tre)
for x in $1/*.fas; do
  o=${x/fas/tre}
  if [[ ! "${DONE[@]}" =~ "${o}" ]]; then
    r=$($HOME/bin/fasttree -nt -gtr -quiet -nopr -gamma $x)
    OUTN+=($o)
    RES+=($r)
    echo $x
  fi
  if [[ ${#RES[@]} == 10 ]]; then
    for (( i = 0 ; i < ${#RES[@]} ; i++ )) do
        echo ${RES[$i]} > ${OUTN[$i]}
    done
    RES=()
    OUTN=()
    echo ${#RES[@]}
  fi
done

for (( i = 0 ; i < ${#RES[@]} ; i++ )) do
        echo ${RES[$i]} > ${OUTN[$i]}
done
