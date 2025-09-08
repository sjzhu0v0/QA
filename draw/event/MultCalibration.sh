# cd /home/szhu/work/alice/analysis/QA/draw/event
# root -l -n -q MultCalibration.cpp


root -l -b << eof
.L MultCalibration.cpp
MultCalibration("/home/szhu/work/alice/analysis/QA/input/event/MultRaw_LHC24pass1_DiElectron.root",
"/home/szhu/work/alice/analysis/QA/output/event/MultCalibration_LHC24pass1_DiElectron.root",
"LHC24pass1_DiElectron" , 35)
eof
