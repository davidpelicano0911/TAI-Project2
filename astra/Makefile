CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra

all: analyse compress decompress

analyse: analyse.cpp
	$(CXX) $(CXXFLAGS) analyse.cpp -o analyse

compress: compress.cpp
	$(CXX) $(CXXFLAGS) compress.cpp -o compress

decompress: decompress.cpp
	$(CXX) $(CXXFLAGS) decompress.cpp -o decompress

test: compress decompress
	@echo "--- Compressing ---"
	@for f in ../data2/*; do \
		b=$$(basename $$f); \
		./compress $$f /tmp/$$b.astr 2>&1 | sed "s/^/  $$b: /"; \
	done
	@echo ""
	@echo "--- Decompressing and verifying ---"
	@ok=0; fail=0; \
	for f in ../data2/*; do \
		b=$$(basename $$f); \
		./decompress /tmp/$$b.astr /tmp/$$b.dec 2>/dev/null; \
		if cmp -s $$f /tmp/$$b.dec; then \
			echo "  $$b: OK"; ok=$$((ok+1)); \
		else \
			echo "  $$b: MISMATCH!"; fail=$$((fail+1)); \
		fi; \
	done; \
	echo ""; echo "Passed: $$ok  Failed: $$fail"

clean:
	rm -f analyse compress decompress
