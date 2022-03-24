.PHONY: all clean

CXX=g++
CXXFLAGS=-std=c++17 -Wall -pedantic

all: strassen

strassen: src/main.cpp
	$(CXX) $(CXXFLAGS) -o $@ -Iinclude $<

clean:
	rm -rf strassen