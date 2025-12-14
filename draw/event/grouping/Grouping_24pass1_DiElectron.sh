# /home/szhu/work/alice/analysis/QA/process/event/Grouping_24pass1_DiElectron.sh

cd /home/szhu/work/alice/analysis/QA/draw/event/grouping
root -l -q -b Grouping_24pass1_DiElectron.cpp

root -l -q -b draw_proportion_main_vs_threshold.cpp