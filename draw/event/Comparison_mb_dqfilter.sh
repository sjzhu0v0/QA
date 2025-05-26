cd /home/szhu/work/alice/analysis/QA/process/event
root -l -q Trigger_Mult.cpp"("\"/home/szhu/work/alice/analysis/QA/input/event/mult_LHC22o_pass4_thin_526641_dqfilter.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC22o_pass4_thin_526641_dqfilter.root\")"

root -l -q Trigger_Mult.cpp"("\"/home/szhu/work/alice/analysis/QA/input/event/mult_LHC22o_pass4_thin_526641_mb.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC22o_pass4_thin_526641_mb.root\")"

cd /home/szhu/work/alice/analysis/QA/draw/event
root -l Comparison_mb_dqfilter.cpp
