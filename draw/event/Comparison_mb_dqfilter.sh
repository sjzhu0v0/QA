# cd /home/szhu/work/alice/analysis/QA/process/event
# root -l -q Trigger_Mult.cpp"("\"/home/szhu/work/alice/analysis/QA/input/event/mult_LHC22o_pass4_thin_526641_dqfilter.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC22o_pass4_thin_526641_dqfilter.root\")"

# root -l -q Trigger_Mult.cpp"("\"/home/szhu/work/alice/analysis/QA/input/event/mult_LHC22o_pass4_thin_526641_mb.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC22o_pass4_thin_526641_mb.root\")"

# cd /home/szhu/work/alice/analysis/QA/draw/event
# root -l Comparison_mb_dqfilter.cpp
# ==============================================================

# cd /home/szhu/work/alice/analysis/QA/process/event
# root -l -q Trigger_Mult.cpp"("\"/home/szhu/work/alice/analysis/QA/input/event/mult_LHC24af_pass1_550367_550632_DiElectron.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC24af_pass1_550367_550632_DiElectron.root\")"

# root -l -q Trigger_Mult.cpp"("\"/home/szhu/work/alice/analysis/QA/input/event/mult_LHC24af_pass1_550367_550632_mb.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC24af_pass1_550367_550632_mb.root\")"

cd /home/szhu/work/alice/analysis/QA/draw/event
root -l Comparison_mb_dqfilter.cpp"("\"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC24af_pass1_550367_550632_mb.root\"", \"/home/szhu/work/alice/analysis/QA/output/event/triggerStudy_LHC24af_pass1_550367_550632_DiElectron.root\", \"_LHC24af_pass1_550367_550632\")"