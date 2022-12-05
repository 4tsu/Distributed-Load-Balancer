MAIN=src/main.cpp
SRCDIR=./src
SRC=$(shell ls $(SRCDIR)/*.cpp)
OBJ=$(SRC:.cpp=.o)
TESTOBJ=$(SRC:.cpp=_test.o)
# 2d simulation
TWODDIR=./src/2d
TWODSRC=$(shell ls $(TWODDIR)/*.cpp)
TWODOBJ=$(TWODSRC:.cpp=.o)
TWODTESTOBJ=$(TWODSRC:.cpp=_test.o)

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



2d: $(TWODOBJ) $(MAIN)
	$(CC) $(OPTIONS) -O3 $^ -o md2d.exe

$(TWODDIR)/%.o: $(TWODDIR)/%.cpp
	$(CC) $(OPTIONS) -O3 -c $< -o $@



# ===test====================================
test.exe: $(TESTOBJ)
	$(CC) $(OPTIONS) $(TESTOPT) $^ -o $@

$(SRCDIR)/%_test.o: $(SRCDIR)/%.cpp
	$(CC) $(OPTIONS) $(TESTOPT) -c $< -o $@

test: test.exe
	-rm *.cdv *.temp 
	-rm energy.dat time_*.dat
	mpirun --oversubscribe -np 4 ./test.exe
	-gnuplot $(VISDIR)/energy_test.plt

dumperr: test.exe
	-rm err.dat
	-rm *.cdv *.temp
	-rm energy.dat time_*.dat
	mpirun --oversubscribe -np 4 ./test.exe 2> err.dat
	-gnuplot $(VISDIR)/energy_test.plt



# 2d simulation
test2d.exe: $(TWODTESTOBJ) $(MAIN)
	$(CC) $(OPTIONS) $(TESTOPT) $^ -o $@

$(TWODDIR)/%_test.o: $(TWODDIR)/%.cpp
	$(CC) $(OPTIONS) $(TESTOPT) -c $< -o $@

test2d: test2d.exe
	-rm *.cdv *.temp 
	-rm energy.dat time_*.dat
	mpirun --oversubscribe -np 4 ./test2d.exe
	-gnuplot $(VISDIR)/energy_test.plt



# ===========================================
dep:
	g++ $(DEPFLAGS) $(SRC) $(TWODSRC) $(OPTIONS) >makefile.depend

run: md.exe
	-rm *.cdv *.temp
	-rm energy.dat time_*.dat
	mpirun -np 4 ./md.exe

fig:
	python3 $(VISDIR)/vis.py

clean:
	rm -f md.exe test.exe md2d.exe test2d.exe
	rm -f $(SRCDIR)/*.o $(TWODDIR)/*.o
	-rm *.cdv *.dat *.temp

# ===========================================