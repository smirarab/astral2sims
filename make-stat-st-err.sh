#!/bin/bash

for x in */*/*.1000; do 
   test -s $x.score && continue; 
   echo -n $x" "; 
   compareTrees.missingBranch ${x%/*}/s_tree.trees $x |tee $x.score; 
done

grep " " */*/*.score |sed -e "s/\// /g" -e "s/main.//g" -e "s/000//" -e "s/\./k /" -e "s/\.1000/ 1000/" -e "s/.score:/ /g"|tee scores.stat
