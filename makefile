SRC=$(shell ls *.cpp)
CC = mpic++
OPTIONS = -std=c++17 -include lib.hpp
OBJ=$(SRC:.cpp=.o)
# for make dep
DEPFLAGS=-MM -MG


.PHONY : all
all: md.exe

md.exe: $(OBJ)
	$(CC) $(OPTIONS) -O3 $(OBJ) -o $@

%.o: %.cpp
	$(CC) $(OPTIONS) -c -O3 $<

test: $(OBJ)
	$(CC) $(OPTIONS) -Wall -Wextra --pedantic-error $(OBJ) -o test.exe
	mpirun -np 4 ./test.exe > e.dat
	-gnuplot energy.plt

run: md.exe
	mpirun -np 4 ./md.exe > e.dat

dep:
	g++ $(DEPFLAGS) $(SRC) $(OPTIONS) >makefile.depend

clean:
	rm -f md.exe *.o