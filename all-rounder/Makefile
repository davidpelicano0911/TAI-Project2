CXX      = g++
CXXFLAGS = -O3 -std=c++17 -pthread -Wall -Wextra
HEADERS  = coder/FseTable.hpp coder/FseEncoder.hpp coder/FseDecoder.hpp \
           model/Predictors.hpp model/CompressorModel.hpp

all: compress decompress

compress: compress.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o compress compress.cpp

decompress: decompress.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o decompress decompress.cpp

clean:
	rm -f compress decompress
