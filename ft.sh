#!/bin/bash

x=$1/all-genes.phylip

o=${x/all-genes.phylip/estimatedgenetre}

$HOME/bin/fasttree -nt -gtr -quiet -nopr -gamma -n 1000 $x > $o

gzip $x

