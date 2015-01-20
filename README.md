Simulation Scripts
===========

Download and unzip [this file](simulation-scripts.zip)

- Use `estimate-parameters.r` to estimate model parameters from values on three different sets of real data (data also included for 1kp, avian, and mammalian datasets)

- Use `simulate.sh` to simulate the species tree and the gene trees (you need to have version 1.0 of [SimPhy](https://code.google.com/p/simphy-project/) installed)

- `simulate.sh` creates control.txt files for each replicate. You need to go to each of these directories and run `indelible` from those directories. Alternatively, use our `run-indelible.sh` script. 

