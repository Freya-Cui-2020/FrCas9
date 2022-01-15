#!/bin/sh
script_dir=PATH 
mkdir Faecalibaculum_rodentium
cd Faecalibaculum_rodentium
mkdir MPJZ01000004.1
cd MPJZ01000004.1
bash $script_dir/pipe_MOD20191111.sh ../../CasPDB_contig/2A/Faecalibaculum_rodentium/ MPJZ01000004.1.fasta `pwd`/MPJZ01000004.1 1> MPJZ01000004.1.out 2> MPJZ01000004.1.err