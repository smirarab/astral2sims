#!/bin/bash

rep=50
genes=1000
b=0.000001
t=500000
sp=200
for ind in 2 5; do  
	./simphy -rs $rep -rl U$genes,$genes -rg 1 -st U$t,$t -si U$ind,$ind -sl U$sp,$sp -sb U$b,$b -p U200000,200000  -hs L1.5,1 -hl L1.2,1 -hg l1.4,1 -u E10000000 -so U1,1 -od 1 -or 0 -v 3  -cs 293745 -o model.$sp-$ind.$t.$b |grep -E "[:-]"| tee log.$sp.$t.$b; 
	rm model.$sp-$ind.$t.$b/*/l_trees.trees
	perl post_stidsim.pl `pwd`/model.$sp-$ind.$t.$b 1 
	for r in `ls -d model.$sp-$ind.$t.$b/*`; do  
		cat $r/g_trees*.trees > $r/truegenetrees; 
		rm  $r/g_trees*.trees; 
	done
done 

