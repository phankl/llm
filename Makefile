#Generic c++ Makefile

CXX = g++
CXXFLAGS = -O0 -g -fopenmp

MAIN = llm

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

$(MAIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(OBJ) $(MAIN)
