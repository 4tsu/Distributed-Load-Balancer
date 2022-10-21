SRC=$(shell ls *.cpp)
CC = mpic++
OPTIONS = -std=c++17 -include lib.hpp
TESTOPT = -Wall -Wextra --pedantic-error
OBJ=$(SRC:.cpp=.o)
# for make dep
DEPFLAGS=-MM -MG


.PHONY : all
all: md.exe

md.exe: $(OBJ)
	$(CC) $(OPTIONS) -O3 $(OBJ) -o $@

%.o: %.cpp
	$(CC) $(OPTIONS) $(TESTOPT) -O3 -c $<

test.exe: $(OBJ)
	$(CC) $(OPTIONS) $(TESTOPT) $(OBJ) -O3 -o $@

test: test.exe
	-rm err.dat
	mpirun --oversubscribe -np 6 ./test.exe > e.dat
	-gnuplot energy.plt

dumperr: test.exe
	-rm err.dat
	mpirun -np 4 ./test.exe > e.dat 2> err.dat

run: md.exe
	mpirun -np 4 ./md.exe > e.dat

dep:
	g++ $(DEPFLAGS) $(SRC) $(OPTIONS) >makefile.depend

clean:
	rm -f md.exe *.o test.exe