CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
TARGET   = benchmark

$(TARGET): benchmark.cpp
	$(CXX) $(CXXFLAGS) benchmark.cpp -o $(TARGET)

run: $(TARGET)
	./$(TARGET) ../data2 -o results.csv

clean:
	rm -f $(TARGET)
