SRCDIR=./src
SRC=$(shell ls $(SRCDIR)/*.cpp)
OBJ=$(SRC:.cpp=.o)

CC = mpic++
OPTIONS = -std=c++17 -include $(SRCDIR)/lib.hpp
TESTOPT = -Wall -Wextra --pedantic-error -Wno-cast-function-type

# for make dep
DEPFLAGS=-MM -MG



all: md.exe
md.exe: $(OBJ)
	$(CC) $(OPTIONS) -O3 $(OBJ) -o $@
	-rm *.cdv

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(OPTIONS) $(TESTOPT) -O3 -c $< -o $@

test.exe: $(OBJ)
	$(CC) $(OPTIONS) $(TESTOPT) $(OBJ) -O3 -o $@

test: test.exe
	-rm *.cdv
	mpirun --oversubscribe -np 4 ./test.exe > e.dat
	-gnuplot energy.plt

dumperr: test.exe
	-rm err.dat
	-rm *.cdv
	mpirun --oversubscribe -np 4 ./test.exe > e.dat 2> err.dat

run: md.exe
	mpirun -np 4 ./md.exe > e.dat

dep:
	g++ $(DEPFLAGS) $(SRC) $(OPTIONS) >makefile.depend

clean:
	rm -f md.exe $(SRCDIR)/*.o test.exe *.cdv