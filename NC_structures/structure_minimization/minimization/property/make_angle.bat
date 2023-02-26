#!/bin/csh -f

set exe = lmp_serial
@ sim = 1
@ st = 3
@ end = $st + $sim - 1
@ end1 = $end + 1
while ($st < $end1)

#echo "$st simulation running"

python pbipb_ang.py ../qsim_cube$st/3_cgmin_max/min.traj

mv avg_pbipb.dat data_pbipb_ang/hist_pbipb.cube$st
mv sur_pbipb.dat data_pbipb_ang/surf_pbipb.cube$st

@ st = $st + 1
end

