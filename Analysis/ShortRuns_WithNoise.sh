#!/bin/bash
#$ -q long
#$ -pe smp 6
#$ -M jnganga@nd.edu
#$ -m abe
module load matlab
matlab -nodisplay -nosplash <ShortRuns_All_Methods_NoisedCntrlStates.m> ShortRuns_WithNoise