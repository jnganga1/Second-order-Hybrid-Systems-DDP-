#!/bin/bash
#$ -q long
#$ -pe smp 6
#$ -M jnganga@nd.edu
#$ -m abe
module load matlab
matlab -nodisplay -nosplash <VaryMethods_TradMod.m> VaryMethods_TradMod