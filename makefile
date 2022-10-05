SRC=$(shell ls *.cpp)
options = -std=c++17 -Wall -Wextra -include lib.hpp

.PHONY : all
all: md.exe

md.exe: $(SRC)
	mpic++ $(options) -O3 $(SRC) -o $@

test: $(SRC)
	mpic++ $(options) --pedantic-error $(SRC) -o test.exe
	./test.exe

run: md.exe
	mpirun -np 4 ./md.exe

clean:
	rm -f md.exe *.o