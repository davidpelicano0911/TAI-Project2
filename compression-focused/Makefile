CXX      = g++
CXXFLAGS = -O3 -std=c++17 -pthread -Wall -Wextra

all: compress decompress range_compress range_decompress

compress: compress.cpp
	$(CXX) $(CXXFLAGS) -o compress compress.cpp

decompress: decompress.cpp
	$(CXX) $(CXXFLAGS) -o decompress decompress.cpp

range_compress: range_compress.cpp
	$(CXX) $(CXXFLAGS) -o range_compress range_compress.cpp

range_decompress: range_decompress.cpp
	$(CXX) $(CXXFLAGS) -o range_decompress range_decompress.cpp

clean:
	rm -f compress decompress range_compress range_decompress
