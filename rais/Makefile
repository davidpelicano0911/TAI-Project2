CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

all: compress decompress coder_compare

compress: compress.cpp
	$(CXX) $(CXXFLAGS) compress.cpp -o compress

decompress: decompress.cpp
	$(CXX) $(CXXFLAGS) decompress.cpp -o decompress

coder_compare: coder_compare.cpp
	$(CXX) $(CXXFLAGS) coder_compare.cpp -o coder_compare

test: compress decompress
	@echo "--- Compressing all data files ---"
	@for f in ../data2/*; do \
		echo -n "$$f -> /tmp/$$(basename $$f).rais  "; \
		./compress $$f /tmp/$$(basename $$f).rais; \
	done
	@echo ""
	@echo "--- Decompressing and verifying ---"
	@ok=0; fail=0; \
	for f in ../data2/*; do \
		b=$$(basename $$f); \
		./decompress /tmp/$$b.rais /tmp/$$b.dec; \
		if cmp -s $$f /tmp/$$b.dec; then \
			echo "$$b: OK"; ok=$$((ok+1)); \
		else \
			echo "$$b: MISMATCH!"; fail=$$((fail+1)); \
		fi; \
	done; \
	echo ""; echo "Passed: $$ok  Failed: $$fail"

clean:
	rm -f compress decompress coder_compare
