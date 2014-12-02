#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -v LD_LIBRARY_PATH
#$ -v PATH
#$ -V

#set -x

H=/home/smirarab/broad_sims/freq_fixed/final_first_1000rep/


test "$SGE_TASK_ID" == "" && SGE_TASK_ID=$1
id=`echo $(printf "%0.4d" $SGE_TASK_ID)`
echo working on $id

method=star

mkdir $H/$id/$method
cd $H/$id/$method

test -s $method.tre && echo output already exists. 
test -s $method.tre && exit 1

/usr/bin/time -a -o $method.time Rscript /home/smirarab/broad_sims/steac_star_njst.r $method 1>$method.out 2>$method.err
