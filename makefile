SRCDIR=./src
SRC=$(shell ls $(SRCDIR)/*.cpp)
OBJ=$(SRC:.cpp=.o)
TESTOBJ=$(SRC:.cpp=_test.o)

CC = mpic++
OPTIONS = -std=c++17 -include $(SRCDIR)/lib.hpp
TESTOPT = -Wall -Wextra --pedantic-error -Wno-cast-function-type -g -O0

# if filesystem is available
OPTIONS += -DFS

# for make dep
DEPFLAGS=-MM -MG

# for visualization
VISDIR=vis
CDVDIR=cdv


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
	mpirun --oversubscribe -np 4 ./test.exe
	-gnuplot $(VISDIR)/energy_test.plt

dumperr: test.exe
	-rm err.dat
	-rm *.cdv
	mpirun --oversubscribe -np 4 ./test.exe 2> err.dat
	-gnuplot $(VISDIR)/energy_test.plt


# ===========================================
dep:
	g++ $(DEPFLAGS) $(SRC) $(OPTIONS) >makefile.depend

run: md.exe
	-rm *.cdv
	mpirun -np 4 ./md.exe

fig:
	python3 $(VISDIR)/vis.py

clean:
	rm -f md.exe $(SRCDIR)/*.o test.exe 
	-rm *.cdv *.dat

# ===========================================