Simulation Scripts
===========

Download and unzip [this file](simulation-scripts.zip)

- Use `R --save < estimate-parameters.r` to estimate model parameters used in the simulation. We estimate these parameters using maximum likelihood (as implemented in R package `sirt`). As input to the ML calculations, we use GTR parameters that we estimated using RAxML for three different sets of real data (see below). The `sirt` package was used to estimate simulation parameters that best capture the distribution across all these three sets of data.
    - [1kp](http://www.pnas.org/content/early/2014/10/28/1323926111): `bases.1p` and `rates.1kp` (518 gene trees based on C12 alignments)
    - [avian](http://www.sciencemag.org/content/346/6215/1320.full): `bases.avaian` and `rates.avian` (14,446 gene trees based on exons, introns, and UCEs)
    - [Song et. al. mammalian dataset](http://www.pnas.org/content/109/37/14942.short): `bases.song` and `rates.song` (424 gene trees)

400 genes are randomly sampled from each dataset. For avian, three sets of 400 genes are sampled from exons, introns, and UCEs, and only restricting to genes between 1000 and 5000bp long. 


- Use `simulate.sh` to simulate the species tree and the gene trees (you need to have version 1.0 of [SimPhy](https://code.google.com/p/simphy-project/) installed, but a MAC binary is included also)

- `simulate.sh` creates `control.txt` files for each replicate. You need to go to each of these directories and run `indelible` from those directories. Alternatively, use our `run-indelible.sh` script. 

