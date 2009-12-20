
PLATFORM := $(shell uname -sri | sed 's/ /-/g')


CXXFLAGS = -g -Wall 

PROG = lpalign
SRCS = main.cpp SW.cpp sequence.cpp alignment.cpp model.cpp lp.cpp shapecontext.cpp 

OBJS = $(SRCS:%.cpp=%.o)
SWIGS = $(SRCS:%.cpp=%_swig.o)
WRAPS = $(SRCS:%.cpp=%_wrap.o)

all: $(PROG) 
$(PROG): main.o $(OBJS)
	g++ -g -o $(PROG) main.o sequence.o  model.o alignment.o lp.o SW.o shapecontext.cpp -L/net/irdev22/export/crawlspace/longbin/install/glpk/lib  -lpthread -lglpk 
SW.o: SW.cpp 
	g++ -g -c -Wall SW.cpp
sequence.o: sequence.cpp
	g++ -g -c -Wall sequence.cpp
model.o: model.h 
	g++ -g -c -Wall model.cpp
alignment.o: alignment.h 
	g++ -g -c -Wall alignment.cpp
lp.o: lp.cpp
	g++ -g -c -Wall -I/net/irdev22/export/crawlspace/longbin/install/glpk/include  lp.cpp
main.o: main.cpp
	g++ -g -Wall -I/net/irdev22/export/crawlspace/longbin/install/glpk/include -c main.cpp
clean:
	rm -rf *~ $(OBJS) $(PROG) ${SWIGS} mss_wrap.cpp mss_wrap.o main.o _mss.so
python: mss_wrap.o sequence_swig.o model_swig.o alignment_swig.o SW_swig.o shapecontext_swig.o
	g++ -shared mss_wrap.o sequence_swig.o model_swig.o alignment_swig.o SW_swig.o shapecontext_swig.o -o _mss.so
shapecontext_swig.o: shapecontext.cpp
	gcc -fPIC -c shapecontext.cpp -o shapecontext_swig.o 
SW_swig.o: SW.cpp
	gcc -fPIC -c SW.cpp -o SW_swig.o 
sequence_swig.o: sequence.cpp
	gcc -fPIC -c sequence.cpp -o sequence_swig.o 
model_swig.o: model.cpp
	gcc -fPIC -c model.cpp -o model_swig.o 
alignment_swig.o: alignment.cpp
	gcc -fPIC -c alignment.cpp -o alignment_swig.o 
mss_wrap.o: mss_wrap.cpp
	gcc -fPIC -c mss_wrap.cpp -o mss_wrap.o -I/usr/include/python2.6
mss_wrap.cpp: mss.i
	swig -c++ -python -o mss_wrap.cpp mss.i 







