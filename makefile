SRC=$(shell ls *.cpp)
CC = ccache mpic++
OPTIONS = -std=c++17 -include lib.hpp

.PHONY : all
all: md.exe

md.exe: $(SRC)
	$(CC) $(OPTIONS) -O3 $(SRC) -o $@

test: $(SRC)
	$(CC) $(OPTIONS) -Wall -Wextra --pedantic-error $(SRC) -o test.exe
	mpirun -np 4 ./test.exe > e.dat
	gnuplot energy.plt

run: md.exe
	mpirun -np 4 ./md.exe > e.dat

clean:
	rm -f md.exe *.o