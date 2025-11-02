#ifndef EventData_h
#define EventData_h

#include "ROOT/RVec.hxx"
#include "iostream"
#include "tuple"
#include "vector"

using namespace std;

struct EventInfo {
  Int_t fMultTPC = 0;
  Int_t fMultTracklets = 0;
  Int_t fMultNTracksPV = 0;
  Float_t fMultFT0C = 0;
  unsigned short fNumContrib = 0;
  Float_t fNumContribCalib = 0; // Calibrated number of contributors
  Float_t fPosX = 0;
  Float_t fPosY = 0;
  Float_t fPosZ = 0;
  Long64_t fSelection = 0;
  Float_t fHadronicRate = 0;

  void Copy(EventInfo other);
  bool isGood() const;
  ClassDefNV(EventInfo, 1); // 添加这一行
};

struct JpsiInfo {
  Int_t fPT_size = 0;
  std::vector<Float_t> fPT;

  Int_t fEta_size = 0;
  std::vector<Float_t> fEta;

  Int_t fPhi_size = 0;
  std::vector<Float_t> fPhi;

  Int_t fMass_size = 0;
  std::vector<Float_t> fMass;

  Int_t fSign_size = 0;
  std::vector<Float_t> fSign;

  template <typename T> void Copy(const T &other) {
    fPT_size = other.fPT_size;
    fPT.assign(other.fPT.begin(), other.fPT.end());
    fEta_size = other.fEta_size;
    fEta.assign(other.fEta.begin(), other.fEta.end());
    fPhi_size = other.fPhi_size;
    fPhi.assign(other.fPhi.begin(), other.fPhi.end());
    fMass_size = other.fMass_size;
    fMass.assign(other.fMass.begin(), other.fMass.end());
    fSign_size = other.fSign_size;
    fSign.assign(other.fSign.begin(), other.fSign.end());
  }

  bool isGood() const;
  ClassDefNV(JpsiInfo, 1); // 添加这一行
};

struct TrackInfo {
  Int_t fPTREF_size = 0;
  std::vector<Float_t> fPTREF;

  Int_t fEtaREF_size = 0;
  std::vector<Float_t> fEtaREF;

  Int_t fPhiREF_size = 0;
  std::vector<Float_t> fPhiREF;

  template <typename T> void Copy(const T &other) {
    fPTREF_size = other.fPTREF_size;
    fPTREF.assign(other.fPTREF.begin(), other.fPTREF.end());
    fEtaREF_size = other.fEtaREF_size;
    fEtaREF.assign(other.fEtaREF.begin(), other.fEtaREF.end());
    fPhiREF_size = other.fPhiREF_size;
    fPhiREF.assign(other.fPhiREF.begin(), other.fPhiREF.end());
  }

  bool isGood() const;
  ClassDefNV(TrackInfo, 1); // 添加这一行
};

struct EventData {
  EventInfo event_info;
  EventInfo event_info2;
  JpsiInfo jpsi_info;
  TrackInfo track_info;

  bool isGood() const;

  ClassDefNV(EventData, 1); // 添加这一行
};

struct EventDataREF {
  EventInfo event_info;
  EventInfo event_info2;
  TrackInfo track_info;
  TrackInfo track_info2;

  ClassDefNV(EventDataREF, 1); // 添加这一行
};

EventData CreateEventData(Int_t fMultTPC, Int_t fMultTracklets,
                          Int_t fMultNTracksPV, Float_t fMultFT0C,
                          unsigned short fNumContrib, double fNumContribCalib,
                          Float_t fPosX, Float_t fPosY, Float_t fPosZ,
                          Long64_t fSelection, Float_t fHadronicRate,
                          const ROOT::VecOps::RVec<Float_t> &fPT,
                          const ROOT::VecOps::RVec<Float_t> &fEta,
                          const ROOT::VecOps::RVec<Float_t> &fPhi,
                          const ROOT::VecOps::RVec<Float_t> &fMass,
                          const ROOT::VecOps::RVec<Float_t> &fSign,
                          const ROOT::VecOps::RVec<Float_t> &fPTREF,
                          const ROOT::VecOps::RVec<Float_t> &fEtaREF,
                          const ROOT::VecOps::RVec<Float_t> &fPhiREF);

#endif // EventData_h