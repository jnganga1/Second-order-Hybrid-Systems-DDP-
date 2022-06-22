#!/bin/bash
#$ -q long
#$ -pe smp 6
#$ -M jnganga@nd.edu
#$ -m abe
module load matlab
matlab -nodisplay -nosplash <ShortRuns_Tens_NoisedCntrlStates.m> ShortRuns_WithNoise