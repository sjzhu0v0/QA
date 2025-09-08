cd /home/szhu/work/alice/analysis/QA/draw/event
# root -l MultCalibrationResult.cpp

root -l -b << eof
.L MultCalibrationResult.cpp
MultCalibrationResult("/home/szhu/work/alice/analysis/QA/input/event/MultCalib_LHC24pass1_DiElectron.root",
"/home/szhu/work/alice/analysis/QA/output/event/MultCalibrationResult_LHC24pass1_DiElectron.root",
"LHC24pass1_DiElectron")
eof