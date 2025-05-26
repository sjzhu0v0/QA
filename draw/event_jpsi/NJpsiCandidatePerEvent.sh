# root -l -b -q /home/szhu/work/alice/analysis/QA/process/event_jpsi/NJpsiCandidatePerEvent.cpp"("\"/home/szhu/data/PairFlow/22pass4_highIR/sample/O2dqflowvecd.root\"", \"/home/szhu/work/alice/analysis/QA/output/event_jpsi/NJpsiCandidatePerEvent.root\")"

# cd /home/szhu/work/alice/analysis/QA/draw/event_jpsi
# root -l NJpsiCandidatePerEvent.cpp"(\"/home/szhu/work/alice/analysis/QA/output/event_jpsi/NJpsiCandidatePerEvent.root\")"

cd /home/szhu/work/alice/analysis/QA/draw/event_jpsi
root -l NJpsiCandidatePerEvent.cpp"(\"/home/szhu/data/DoubleJpsi/22pass4_highIR/doubleJpsi.root\")"