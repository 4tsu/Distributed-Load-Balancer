SRC=$(shell ls *.cpp)
OBJ=$(SRC:.cpp=.o)
options=-std=c++17 -Wall -Wextra

all: md.exe

md.exe: $(OBJ)
	mpic++ $(options) -O3 $(OBJ) -o $@

test: $(OBJ)
	mpic++ $(options) --pedantic-error $(OBJ) -o test.exe
	./test.exe

run: md.exe
	mpirun -np 4 ./md.exe

clean:
	rm -f md.exe *.o