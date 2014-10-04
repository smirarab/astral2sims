#!/bin/bash

rep=50
genes=1000

sp=200

for b in 0.000001 0.0000001; do
	for t in 500000 2000000 10000000; do  
		~/workspace/simphy-project/src/simphy -rs $rep -rl U$genes,$genes -rg 1 -st U$t,$t -si U1,1 -sl U$sp,$sp -sb U$b,$b -p U200000,200000  -hs L1.5,1 -hl L1.2,1 -hg l1.4,1 -u E10000000 -so U1,1 -od 1 -or 0 -v 3  -cs 293745 -o main.$t.$b |grep -E "[:-]"| tee log.$t.$b; 
		rm main.$t.$b/*/l_trees.trees
		for r in `ls -d main.$t.$b/*`; do 
			sed -i "" -e "s/_0_0//g" $r/g_trees*.trees; 
		done
		perl post_stidsim.pl `pwd`/main.$t.$b 1 
		#for r in `ls -d main.$t.$b/*`; do  
		#	cat $r/g_trees*.trees > $r/allgenetrees; 
		#	rm  $r/g_trees*.trees; 
		#done
	done; 
done


# find . -name control.txt -execdir sh -c 'condor_run indelible &' \;
