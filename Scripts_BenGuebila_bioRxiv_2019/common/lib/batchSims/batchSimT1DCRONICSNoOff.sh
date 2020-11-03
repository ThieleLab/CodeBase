#!/bin/bash
matlab='/usr/local/MATLAB/R2016b/bin/./matlab'

runList='IVGTT IVITT Noinf SCIB SCII WBLiquid WBSolid' 
INPUTFILE=CRONICSnoOff
OUTPUTFILE=CRONICSnoOff.out
for run in $runList; do(
	foo="../TrialSimulation/$run/ODE_T1D"
	cd "$foo"
	$matlab -nodesktop -nosplash -nodisplay -r $INPUTFILE -logfile $OUTPUTFILE 
) &
done
