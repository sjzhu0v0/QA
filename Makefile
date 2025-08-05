DIR_BASE=/lustre/alice/users/szhu/work/Analysis/PairFlow
PATH_INCLUDE=$(DIR_BASE)/include
FLAGS_INCLUDE=-I$(DIR_BASE)/include -I$(DIR_BASE)/macro
FLAGS_ROOT=$(shell root-config --cflags --libs)
FLAGS_MINUIT=-lMinuit

all: \
	macro/track/TrackInfoMC.exe \
	macro/event/MultRaw.exe \
	macro/event/MultCalib.exe \
	macro/event/MultPileupCut.exe \
	macro/event/MultQA_AllCut.exe \
	macro/event/MultREFRaw.exe \
	macro/jpsi/JpsiQA.exe \
	macro/event_jpsi/NJpsiCandidatePerEvent.exe \
	macro/event_jpsi/EventMixingJpsiAsso.exe \
	macro/event_jpsi/EventMixingJpsiAsso_v2.exe \
	macro/event_jpsi/MixEventReading.exe \
	macro/event_jpsi/JpsiAsso.exe

macro/track/TrackInfoMC.exe: macro/track/TrackInfoMC.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event/MultRaw.exe: macro/event/MultRaw.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event/MultCalib.exe: macro/event/MultCalib.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event/MultPileupCut.exe: macro/event/MultPileupCut.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event_jpsi/NJpsiCandidatePerEvent.exe: macro/event_jpsi/NJpsiCandidatePerEvent.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event/MultQA_AllCut.exe: macro/event/MultQA_AllCut.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event/MultREFRaw.exe: macro/event/MultREFRaw.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT)

macro/event_jpsi/JpsiAsso.exe: macro/event_jpsi/JpsiAsso.cpp
	g++ -o $@ $^ $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT) -I./ -L./opt -lMRootDict

macro/jpsi/JpsiQA.exe: macro/jpsi/JpsiQA.cpp opt/libMRootDict.so opt/libMRootDict.so
	g++ -o $@ macro/jpsi/JpsiQA.cpp $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT) -I./ -L./opt -lMRootDict

macro/event_jpsi/EventMixingJpsiAsso.exe: macro/event_jpsi/EventMixingJpsiAsso.cpp opt/libMRootDict.so
	g++ -o $@ macro/event_jpsi/EventMixingJpsiAsso.cpp $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT) -I./ -L./opt -lMRootDict

macro/event_jpsi/EventMixingJpsiAsso_v2.exe: macro/event_jpsi/EventMixingJpsiAsso_v2.cpp opt/libMRootDict.so
	g++ -o $@ macro/event_jpsi/EventMixingJpsiAsso_v2.cpp $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT) -I./ -L./opt -lMRootDict

macro/event_jpsi/MixEventReading.exe: macro/event_jpsi/MixEventReading.cpp opt/libMRootDict.so
	g++ -o $@ macro/event_jpsi/MixEventReading.cpp $(FLAGS_INCLUDE) $(FLAGS_ROOT) $(FLAGS_MINUIT) -I./ -L./opt -lMRootDict

opt/MRootDict.cxx: opt/EventData.h opt/LinkDef.h
	rootcint -f $@ -c opt/EventData.h opt/LinkDef.h $(FLAGS_INCLUDE)

opt/libMRootDict.so: opt/MRootDict.cxx opt/EventData.cxx
	g++ -o $@ $^ -I./ `root-config --cflags --libs` -shared -fPIC