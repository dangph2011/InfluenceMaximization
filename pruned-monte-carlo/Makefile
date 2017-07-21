CXX = g++
CXXFLAGS = -g -Wall -Wextra -O3 -std=c++11
CXXFLAGS_DEBUG = -ggdb -std=c++11
#INCLUDE = /home/dangph/usr/include

all: benchmark

benchmark: samples/benchmark.cpp src/pmc.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

benchmark_hub: samples/benchmark.cpp src/pmc_hub.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

benchmark_hub_origin: samples/benchmark.cpp src/pmc_origin.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

benchmark_estimate: samples/benchmark.cpp src/pmc_estimate.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

debug: samples/benchmark.cpp src/pmc.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

debug_estimate: samples/benchmark.cpp src/pmc_estimate.cpp src/evaluate.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ 

.PHONY: clean

clean:
	rm -f benchmark benchmark_hub
