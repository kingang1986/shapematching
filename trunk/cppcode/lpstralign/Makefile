
PLATFORM := $(shell uname -sri | sed 's/ /-/g')


CXXFLAGS = -g -Wall 

PROG = lpalign
SRCS = main.cpp SW.cpp sequence.cpp Data.cpp lp.cpp 

OBJS = $(SRCS:%.cpp=%.o)
SWIGS = $(SRCS:%.cpp=%_swig.o)
WRAPS = $(SRCS:%.cpp=%_wrap.o)

all: $(PROG) 
$(PROG): main.o $(OBJS)
	g++ -g -o $(PROG) main.o sequence.o Data.o lp.o SW.o -L/net/irdev22/export/crawlspace/longbin/install/glpk/lib  -lpthread -lglpk 
SW.o: SW.cpp Data.h
	g++ -g -c -Wall SW.cpp
sequence.o: sequence.cpp
	g++ -g -c -Wall sequence.cpp
Data.o: Data.cpp sequence.h
	g++ -g -c -Wall Data.cpp
lp.o: lp.cpp
	g++ -g -c -Wall -I/net/irdev22/export/crawlspace/longbin/install/glpk/include  lp.cpp
main.o: main.cpp
	g++ -g -Wall -I/net/irdev22/export/crawlspace/longbin/install/glpk/include -c main.cpp
clean:
	rm -rf *~ $(OBJS) $(PROG) ${SWIGS} ${WRAPS} main.o _sequence.so
python: sequence_wrap.o sequence_swig.o   
	g++ -shared sequence_wrap.o sequence_swig.o -o _sequence.so
sequence_wrap.o: sequence_wrap.cpp
	gcc -fPIC -c sequence_wrap.cpp -o sequence_wrap.o -I/usr/include/python2.6
sequence_swig.o: sequence.cpp
	gcc -fPIC -c sequence.cpp -o sequence_swig.o 
sequence_wrap.cpp: sequence.i
	swig -c++ -python -o sequence_wrap.cpp sequence.i







