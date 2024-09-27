cd test/schan
../../../ZTOPschan.x < ls_nc.in > ls_nc.log
../../../ZTOPschan.x < ls_tc.in > ls_tc.log
../../../ZTOPschan.x < s_nc.in > s_nc.log
../../../ZTOPschan.x < s_dc.in > s_dc.log
../../../ZTOPschan.x < s_jc.in > s_jc.log
../../../ZTOPschan.x < s_tc.in > s_tc.log
cd ../tchan
../../../ZTOPtchan.x < lt_nc.in > lt_nc.log
../../../ZTOPtchan.x < lt_tc.in > lt_tc.log
../../../ZTOPtchan.x < t_nc.in > t_nc.log
../../../ZTOPtchan.x < t_dc.in > t_dc.log
../../../ZTOPtchan.x < t_jc.in > t_jc.log
../../../ZTOPtchan.x < t_tc.in > t_tc.log
cd ../../
echo 'Test done.  Compare distributions in test/schan & test/tchan to schan & tchan'
