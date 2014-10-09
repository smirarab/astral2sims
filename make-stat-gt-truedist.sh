#!/bin/bash

p=`pwd`;

for x in */0*/; do cd $p/$x && pwd && compareTrees.missingBrach s_tree.trees truegenetrees done |tee gt-truedist.stat

grep " " */*/gt-truedist.stat|sed -e "s/.gt.*:/ /g" -e "s/main.//g" -e "s/\./ /" -e "s/\// /g" -e "s/000//" -e "s/ /k /"|tee gt-err.stat
