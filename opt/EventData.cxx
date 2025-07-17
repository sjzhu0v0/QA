#include "EventData.h"

void EventInfo::Copy(EventInfo other) {
  fMultTPC = other.fMultTPC;
  fMultTracklets = other.fMultTracklets;
  fMultNTracksPV = other.fMultNTracksPV;
  fMultFT0C = other.fMultFT0C;
  fNumContrib = other.fNumContrib;
  fNumContribCalib = other.fNumContribCalib; // Calibrated number of contributors
  fPosX = other.fPosX;
  fPosY = other.fPosY;
  fPosZ = other.fPosZ;
  fSelection = other.fSelection;
  fHadronicRate = other.fHadronicRate;
}

EventData CreateEventData(Int_t fMultTPC, Int_t fMultTracklets,
                          Int_t fMultNTracksPV, Float_t fMultFT0C,
                          unsigned short fNumContrib, Float_t fNumContribCalib,
                          Float_t fPosX, Float_t fPosY, Float_t fPosZ,
                          unsigned long long fSelection, Float_t fHadronicRate,
                          const ROOT::VecOps::RVec<Float_t> &fPT,
                          const ROOT::VecOps::RVec<Float_t> &fEta,
                          const ROOT::VecOps::RVec<Float_t> &fPhi,
                          const ROOT::VecOps::RVec<Float_t> &fMass,
                          const ROOT::VecOps::RVec<Float_t> &fSign,
                          const ROOT::VecOps::RVec<Float_t> &fPTREF,
                          const ROOT::VecOps::RVec<Float_t> &fEtaREF,
                          const ROOT::VecOps::RVec<Float_t> &fPhiREF) {
  EventData event;

  // 填充固定大小的变量
  event.event_info.fMultTPC = fMultTPC;
  event.event_info.fMultTracklets = fMultTracklets;
  event.event_info.fMultNTracksPV = fMultNTracksPV;
  event.event_info.fMultFT0C = fMultFT0C;
  event.event_info.fNumContrib = fNumContrib;
  event.event_info.fNumContribCalib = fNumContribCalib;
  event.event_info.fPosX = fPosX;
  event.event_info.fPosY = fPosY;
  event.event_info.fPosZ = fPosZ;
  event.event_info.fSelection = fSelection;
  event.event_info.fHadronicRate = fHadronicRate;

  // 填充动态数组
  event.jpsi_info.fPT_size = fPT.size();
  event.jpsi_info.fPT.assign(fPT.begin(), fPT.end());

  event.jpsi_info.fEta_size = fEta.size();
  event.jpsi_info.fEta.assign(fEta.begin(), fEta.end());

  event.jpsi_info.fPhi_size = fPhi.size();
  event.jpsi_info.fPhi.assign(fPhi.begin(), fPhi.end());

  event.jpsi_info.fMass_size = fMass.size();
  event.jpsi_info.fMass.assign(fMass.begin(), fMass.end());

  event.jpsi_info.fSign_size = fSign.size();
  event.jpsi_info.fSign.assign(fSign.begin(), fSign.end());

  event.track_info.fPTREF_size = fPTREF.size();
  event.track_info.fPTREF.assign(fPTREF.begin(), fPTREF.end());

  event.track_info.fEtaREF_size = fEtaREF.size();
  event.track_info.fEtaREF.assign(fEtaREF.begin(), fEtaREF.end());

  event.track_info.fPhiREF_size = fPhiREF.size();
  event.track_info.fPhiREF.assign(fPhiREF.begin(), fPhiREF.end());

  return event;
};