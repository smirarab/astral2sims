#!/bin/bash

test `cat gt-err.stat|wc -l` ==  `cat allgenetrees|wc -l` && exit 1;

gsplit --numeric-suffixes=1 -l1 -a4 allgenetrees true.gt.

for x in *.tre; do 
	compareTrees.missingBranch true.gt.${x/.tre/} $x; 
done |tee gt-err.stat

rm true.gt.????
