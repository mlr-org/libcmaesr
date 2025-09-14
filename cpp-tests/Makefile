# Incremental build Makefile for libcmaes tests

# Absolute project root, derived from this Makefile's location
TOP := $(abspath ..)

# libcmaes paths
LIBCMAES_SRC_DIR := $(TOP)/src/libcmaes
LIBCMAES_BUILD_DIR := $(LIBCMAES_SRC_DIR)/build
LIBCMAES_LIB := $(LIBCMAES_BUILD_DIR)/src/libcmaes.a

# Eigen include dir
EIGEN_DIR := /usr/include/eigen3

# Compiler and flags
CXX ?= g++
CXXFLAGS ?= -std=c++11  \
  -I$(EIGEN_DIR) \
  -I$(LIBCMAES_SRC_DIR)/include \
  -I$(LIBCMAES_BUILD_DIR)/include

LDFLAGS ?= -L$(LIBCMAES_BUILD_DIR)/src -Wl,-Bstatic -lcmaes -Wl,-Bdynamic

# CMake configuration flags
CMAKE_FLAGS := \
  -DEIGEN3_INCLUDE_DIR=$(EIGEN_DIR) \
  -DLIBCMAES_USE_OPENMP=OFF \
  -DLIBCMAES_BUILD_SHARED_LIBS=OFF \
  -DLIBCMAES_ENABLE_SURROG=OFF

# Test binaries
TESTS := test_ipop_max_fevals test_bipop_max_fevals test_fevals_mismatch

.PHONY: all ipop bipop libcmaes clean clean-lib

all: $(TESTS)

ipop: test_ipop_max_fevals

bipop: test_bipop_max_fevals

# Ensure the build directory is configured, then build the static library target
$(LIBCMAES_LIB):
	@mkdir -p $(LIBCMAES_BUILD_DIR)
	cmake -S $(LIBCMAES_SRC_DIR) -B $(LIBCMAES_BUILD_DIR) $(CMAKE_FLAGS)
	cmake --build $(LIBCMAES_BUILD_DIR) --target cmaes -j

# Link tests against the freshly built libcmaes

test_ipop_max_fevals: test_ipop_max_fevals.cpp $(LIBCMAES_LIB)
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) -o $@
	./test_ipop_max_fevals

test_bipop_max_fevals: test_bipop_max_fevals.cpp $(LIBCMAES_LIB)
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) -o $@
	./test_bipop_max_fevals

test_fevals_mismatch: test_fevals_mismatch.cpp $(LIBCMAES_LIB)
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) -o $@
	./test_fevals_mismatch

# Cleaning
clean:
	rm -f $(TESTS)

clean-lib:
	rm -rf $(LIBCMAES_BUILD_DIR)
