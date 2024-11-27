#!/bin/csh
#

cd $TMPDIR

module load molpro/2020.1/openmp

set point = $1
set endPoint = $2
echo $point 
while ($point < $endPoint)
    set line = `awk 'NR=='"$point" /home/zcaposm/Scratch/Methanol/MOLPRO/Threshold/2D/Grid2D.txt`

    if ("$line" != "") then
        csh -f /home/zcaposm/Scratch/Methanol/MOLPRO/Threshold/2D/GenerateMolproScript2D.csh $line
    endif
@ point++
end