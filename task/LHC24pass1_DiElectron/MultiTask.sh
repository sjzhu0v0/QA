path_macro="/lustre/alice/users/szhu/work/Analysis/QA"
path_output="/lustre/alice/users/szhu/job/QA/LHC24pass1_DiElectron_Group"
path_input="/lustre/alice/users/szhu/job/JpsiFlowPair24/output/"
path_calib="/lustre/alice/users/szhu/work/Analysis/InfoRun/MultCalib/MultCalibration_LHC24pass1_DiElectron.root:fNumContribfPosZRun_calib_"
tag1=${1}

${path_macro}/macro/jpsi/JpsiQA.exe ${path_input}/${tag1}/tree.root ${path_output}/JpsiQA_${tag1}.root ${tag1} ${path_calib}
${path_macro}/macro/event/MultREFRaw.exe ${path_input}/${tag1}/tree.root ${path_output}/MultREFRaw_${tag1}.root ${tag1} ${path_calib}
${path_macro}/macro/event_jpsi/JpsiAsso.exe ${path_input}/${tag1}/tree.root ${path_output}/JpsiAsso_${tag1}.root ${tag1} 1 ${path_calib}
${path_macro}/macro/event_jpsi/EventMixingJpsiAsso.exe ${path_input}/${tag1}/tree.root ${path_output}/EventMixingJpsiAsso_${tag1}.root ${path_output}/EventMixingJpsiAsso_tree_${tag1}.root ${tag1} ${path_calib}
${path_macro}/macro/event_jpsi/MixEventReading.exe ${path_output}/EventMixingJpsiAsso_tree_${tag1}.root ${path_output}/MixEventReading_${tag1}.root