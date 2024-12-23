#!/bin/csh
#

cd $TMPDIR

module load molpro/2020.1/openmp

set point = $1
set endPoint = $2
echo $point 
while ($point < $endPoint)
    set line = `awk 'NR=='"$point" /home/zcaposm/Scratch/Methanol/MOLPRO/C3V/3D/Grid3D.txt`

    if ("$line" != "") then
        csh -f /home/zcaposm/Scratch/Methanol/MOLPRO/C3V/3D/GenerateMolproScript3D.csh $line
    endif
@ point++
end