CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
TARGET   = benchmark
DATA_DIR ?= ../data2
TMP_DIR  ?= /tmp/tai_bench
CSV      ?= results.csv
PLOT     ?= results_plot.png
BENCHMARK_TOOLS = gzip bzip2 xz zstd

.PHONY: all run csv benchmark-run current check-tools plot clean

all: $(TARGET)

$(TARGET): benchmark.cpp
	$(CXX) $(CXXFLAGS) benchmark.cpp -o $(TARGET)

current:
	$(MAKE) -C ../hais
	$(MAKE) -C ../hais2
	$(MAKE) -C ../helix
	$(MAKE) -C ../rais
	$(MAKE) -C ../astra
	$(MAKE) -C ../prism
	$(MAKE) -C ../prism_block
	$(MAKE) -C ../apex

check-tools:
	@missing=0; \
	for tool in $(BENCHMARK_TOOLS); do \
		if ! command -v $$tool >/dev/null 2>&1; then \
			echo "Missing benchmark tool: $$tool"; \
			missing=1; \
		fi; \
	done; \
	exit $$missing

run csv benchmark-run: $(TARGET) current check-tools
	./$(TARGET) $(DATA_DIR) $(TMP_DIR) -o $(CSV)

plot: $(CSV)
	python3 plot.py $(CSV) $(PLOT)

clean:
	rm -f $(TARGET) $(CSV) $(PLOT)
