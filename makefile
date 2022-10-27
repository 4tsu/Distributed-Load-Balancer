SRCDIR=./src
SRC=$(shell ls $(SRCDIR)/*.cpp)
OBJ=$(SRC:.cpp=.o)
TESTOBJ=$(SRC:.cpp=_test.o)

CC = mpic++
OPTIONS = -std=c++17 -include $(SRCDIR)/lib.hpp
TESTOPT = -Wall -Wextra --pedantic-error -Wno-cast-function-type -g -O0

# for make dep
DEPFLAGS=-MM -MG



# ===release=================================
all: md.exe

md.exe: $(OBJ)
	$(CC) $(OPTIONS) -O3 $^ -o $@

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(OPTIONS) -O3 -c $< -o $@


# ===test====================================
test.exe: $(TESTOBJ)
	$(CC) $(OPTIONS) $(TESTOPT) $^ -o $@

$(SRCDIR)/%_test.o: $(SRCDIR)/%.cpp
	$(CC) $(OPTIONS) $(TESTOPT) -c $< -o $@

test: test.exe
	-rm *.cdv
	mpirun --oversubscribe -np 4 ./test.exe > e.dat
	-gnuplot energy.plt

dumperr: test.exe
	-rm err.dat
	-rm *.cdv
	mpirun --oversubscribe -np 4 ./test.exe > e.dat 2> err.dat
	-gnuplot energy.plt


# ===========================================
run: md.exe
	-rm *.cdv
	mpirun -np 4 ./md.exe > e.dat

dep:
	g++ $(DEPFLAGS) $(SRC) $(OPTIONS) >makefile.depend

clean:
	rm -f md.exe $(SRCDIR)/*.o test.exe *.cdv

# ===========================================