#!/bin/bash

p=`pwd`;

for x in model.*/$1/; do cd $p/$x && pwd && compareTrees.missingBranch s_tree.trees truegenetrees|tee gt-truedist.stat; done 

grep " " */*/gt-truedist.stat|sed -e "s/.gt.*:/ /g" -e "s/model.//g" -e "s/\./ /" -e "s/\// /g" -e "s/000//" -e "s/ /k /"|tee gt-truedist.stat

