cd /home/szhu/work/alice/analysis/QA/process/event

root -l -b << EOF
.L Grouping.cpp
CreateMatrixRun24()
EOF