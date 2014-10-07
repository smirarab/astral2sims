#!/bin/bash

cd $1

TGT=truegenetrees
EST=estimatedgenetrees

test `cat gt-err.stat|wc -l` ==  `cat $TGT|wc -l` && exit 1;

if [ -n "`which gsplit`" ]; then
 gsplit --numeric-suffixes=1 -l1 -a4 $TGT true.gt.
else
 split -d -l1 -a4 $TGT true.gt.
 split -d -l1 -a4 $EST estimated.
fi

for x in estimated.*; do 
	fn=$( compareTrees.missingBranch true.gt.${x/estimated./} $x ); 
	fp=$( compareTrees.missingBranch $x true.gt.${x/estimated./} ); 
        echo ${x/estimated./} $fn $fp
done |tee gt-err.stat

rm estimated.????
