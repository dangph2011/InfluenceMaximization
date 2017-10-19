CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11
CXXFLAGS_DEBUG = -ggdb -std=c++11
#INCLUDE = /home/dangph/usr/include

all: benchmark

benchmark: samples/benchmark.cpp src/pmc.cpp SIEA/GLib.hpp SIEA/mappedHeap.hpp SIEA/HeapData.hpp SIEA/option.cpp SIEA/rwgraph.cpp
	$(CXX) -fopenmp $(CXXFLAGS) SIEA/sfmt/SFMT.c $^ -o $@

benchmark_hub: samples/benchmark.cpp src/pmc_hub.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

benchmark_origin: samples/benchmark.cpp src/pmc_origin.cpp SIEA/GLib.hpp SIEA/mappedHeap.hpp SIEA/HeapData.hpp SIEA/option.cpp SIEA/rwgraph.cpp
	$(CXX) -fopenmp $(CXXFLAGS) SIEA/sfmt/SFMT.c $^ -o $@

benchmark_estimate: samples/benchmark.cpp src/pmc_estimate.cpp SIEA/GLib.hpp SIEA/mappedHeap.hpp SIEA/HeapData.hpp SIEA/option.cpp SIEA/rwgraph.cpp
	$(CXX) -fopenmp $(CXXFLAGS) SIEA/sfmt/SFMT.c $^ -o $@

bin:
	g++ SIEA/el2bin_forward.cpp -o el2bin_forward $(CXXFLAGS)

debug: samples/benchmark.cpp src/pmc.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

debug_estimate: samples/benchmark.cpp src/pmc_estimate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f benchmark benchmark_hub el2bin_forward outinfest
