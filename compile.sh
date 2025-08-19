DIR_BASE=/lustre/alice/users/szhu/work/Analysis/PairFlow
PATH_INCLUDE=${DIR_BASE}/include
FLAGS_INCLUDE="-I${DIR_BASE}/include -I${DIR_BASE}/macro"
FLAGS_ROOT=$(root-config --cflags --libs)
FLAGS_MINUIT=-lMinuit

g++ -o $1 $2 ${FLAGS_INCLUDE} ${FLAGS_ROOT} ${FLAGS_MINUIT}