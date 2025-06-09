DIR_BASE=/lustre/alice/users/szhu/work/Analysis/PairFlow
FLAGS_INCLUDE=-I$(DIR_BASE)/include -I$(DIR_BASE)/macro
FLAGS_ROOT=$(shell root-config --cflags --libs)
FLAGS_MINUIT=-lMinuit

all: \
	macro/event/MultRaw.exe \
	macro/event/MultCalib.exe \
	macro/event/MultPileupCut.exe \
	macro/event_jpsi/NJpsiCandidatePerEvent.exe \
	macro/event/MultQA_AllCut.exe \
	macro/event/MultREFRaw.exe

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