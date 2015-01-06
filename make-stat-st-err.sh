#!/bin/bash

for x in `ls */*/*.genes1000 */*/*.genes200 */*/*.genes50 */*/*.genes*-*ls |grep -v RAxML`; do 
   test -s $x.score && continue; 
   echo -n $x" "; 
   test -s $x || continue; 
   echo `compareTrees.missingBranch ${x%/*}/s_tree.trees $x` `compareTrees.missingBranch $x ${x%/*}/s_tree.trees` |tee $x.score; 
done

grep " " */*/*.score |sed -e "s/\// /g" -e "s/model.//g" -e"s/\./ /" -e "s/000\./K /" -e "s/000K /M /" -e "s/\.genes\(.*\)-\(.*\)ls/ \1 \2/" -e "s/\.genes\(.[0-9]*\)\./ \1 -1 /" -e "s/.score:/ /g"| sed -e "s/.halfresolved//g" > scores.stat
