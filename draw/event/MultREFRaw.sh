cd /home/szhu/work/alice/analysis/QA/draw/event
# root -l MultREFRaw.cpp

root -l -b << eof
.L MultREFRaw.cpp
MultREFRaw("/home/szhu/work/alice/analysis/QA/input/event/MultREFRaw_LHC24_pass1_DiElectron.root", "/home/szhu/work/alice/analysis/QA/output/event/MultREFRaw_LHC24pass1_DiElectron.root", "LHC24pass1_DiElectron")
eof