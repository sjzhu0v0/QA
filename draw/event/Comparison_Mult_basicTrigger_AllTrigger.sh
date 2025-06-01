#  auto rdf_fullTrigger =
#       rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
#           .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
#           .Filter("isntTimeFrameBorder", "no Time Frame border")
#           .Filter("isntSameBunchPileup", "no same bunch pileup")
#           .Filter("isntSelfDefinedPileup", "no self defined pileup");
#   auto rdf_basicTrigger =
#       rdf_witTrigger.Filter("isTriggerTVX", "is Trigger TVX")
#           .Filter("isntITSROFrameBorder", "no ITS RO Frame border")
#           .Filter("isntTimeFrameBorder", "no Time Frame border")
#           .Filter("isntSameBunchPileup", "no same bunch pileup");
cd /home/szhu/work/alice/analysis/QA/draw/event
root -l Comparison_Mult_basicTrigger_AllTrigger.cpp