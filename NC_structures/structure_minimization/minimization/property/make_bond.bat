#!/bin/csh -f

set exe = lmp_serial
@ sim = 1
@ st = 2
@ end = $st + $sim - 1
@ end1 = $end + 1
while ($st < $end1)

#echo "$st simulation running"

python pbi_dist.py ../qsim_cube$st/3_cgmin_max/min.traj

mv avg_pbi.dat data_pbi_dist/hist_pbi.cube$st
mv sur_pbi.dat data_pbi_dist/surf_pbi.cube$st

@ st = $st + 1
end

