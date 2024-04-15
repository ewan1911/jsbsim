rm Zzz_Alt.txt Zzz_aR.txt Zzz_direction.txt Zzz_East.txt Zzz_North.txt Zzz_Roll.txt
rm Zzz_rollError.txt Zzz_rudder.txt Zzz_t.txt Zzz_Time.txt Zzz_timeTurn.txt Zzz_yaw.txt
cd build
make install
cd ..
JSBSim --script=scripts/sgs_test