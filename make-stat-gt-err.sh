#!/bin/bash

p=`pwd`;

for x in model.*/$1/; do cd $p/$x && pwd && /projects/sate4/smirarab/indelible//astral2sims/gt.err.sh `pwd`; done

cd $p
grep " " */*/gt-err.stat|sed -e "s/.gt.*:/ /g" -e "s/model.//g" -e "s/\./ /" -e "s/\./ /" -e "s/\// /g" -e "s/0000 /0k /"|tee gt-err.stat
