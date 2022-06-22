#!/bin/bash
#$ -q long
#$ -pe smp 6
#$ -M jnganga@nd.edu
#$ -m abe
module load matlab
matlab -nodisplay -nosplash <FullRuns_All_Methods_NoisedCntrlStates.m> FullRuns_WithNoise