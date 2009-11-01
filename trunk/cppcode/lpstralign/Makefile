
PLATFORM := $(shell uname -sri | sed 's/ /-/g')


CXXFLAGS = -g -Wall 

PROG = lpalign
SRCS = main.cpp SW.cpp Data.cpp lp.cpp 

OBJS = $(SRCS:%.cpp=%.o)

all: $(PROG) 
$(PROG): main.o $(OBJS)
	g++  -o $(PROG) main.o Data.o lp.o SW.o -L/net/irdev22/export/crawlspace/longbin/install/glpk/lib  -lpthread -lglpk 
SW.o: SW.cpp Data.h
	g++ -c -Wall SW.cpp
Data.o: Data.cpp
	g++ -c -Wall Data.cpp
lp.o: lp.cpp
	g++ -c -Wall -I/net/irdev22/export/crawlspace/longbin/install/glpk/include  lp.cpp
main.o: main.cpp
	g++  -Wall -I/net/irdev22/export/crawlspace/longbin/install/glpk/include -c main.cpp
clean:
	rm -rf *~ $(OBJS) $(PROG) main.o 


