CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra

all: compress decompress

compress: compress.cpp
	$(CXX) $(CXXFLAGS) -o compress compress.cpp

decompress: decompress.cpp
	$(CXX) $(CXXFLAGS) -o decompress decompress.cpp

clean:
	rm -f compress decompress
