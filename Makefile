UNAME_S := $(shell uname -s)

CXX = mpicxx
CXXFLAGS ?= -std=c++20 -fopenmp 
CPPFLAGS := $(CPPFLAGS) -I./include -I./include/muparser -I./include/parameters -I./include/solver -I./include/matrix
SRCS = $(wildcard ./src/*/*.cpp) src/main.cpp
SRCS_TEST = $(wildcard ./src/*/*.cpp) test/main.cpp test/test.cpp 

DOXYGEN ?= doxygen
DOXYFILE ?= Doxyfile

EXEC ?= executable
DOCS ?= documentation
RM ?= rm

MUPARSER_URL = https://github.com/beltoforion/muparser
MUPARSER_DIR = ./muparser
JSON_URL = https://github.com/nlohmann/json
JSON_DIR = ./json

.PHONY = all debug test clean doc libs symlinks allclean

ifeq ($(UNAME_S),Darwin)
	LLVM_PREFIX := $(shell brew --prefix llvm)
    LIBOMP_PREFIX := $(shell brew --prefix libomp)
    OPENMPI_PREFIX := $(shell brew --prefix open-mpi)
	HB_CLANG_CXX := $(LLVM_PREFIX)/bin/clang++
    HB_CLANG_CC := $(LLVM_PREFIX)/bin/clang
    HB_MPICXX := $(OPENMPI_PREFIX)/bin/mpicxx
	CXX := OMPI_CXX=$(HB_CLANG_CXX) $(HB_MPICXX)
	CMAKE_CC := $(HB_CLANG_CC)
    CMAKE_CXX := $(HB_CLANG_CXX)
	CPPFLAGS += -I$(LIBOMP_PREFIX)/include -I$(LLVM_PREFIX)/include
    LDFLAGS += -L$(LIBOMP_PREFIX)/lib -lomp -L$(LLVM_PREFIX)/lib
    LIBS = libmuparser.dylib libmuparser.2.dylib ./json/single_include/nlohmann/json.hpp
else ifeq ($(UNAME_S),Linux)
	CXX := mpicxx
    CMAKE_CC := $(shell which gcc)
    CMAKE_CXX := $(shell which g++)
    LDFLAGS += -fopenmp
    LIBS = libmuparser.so libmuparser.so.2 ./json/single_include/nlohmann/json.hpp
else ifeq ($(UNAME_S),MINGW)
	CXX := mpicxx
    CMAKE_CC := $(shell which gcc)
    CMAKE_CXX := $(shell which g++)
    LDFLAGS += -fopenmp
    LIBS = muparser.dll ./json/single_include/nlohmann/json.hpp
else ifeq ($(UNAME_S),Windows)
	CXX := mpicxx
    CMAKE_CC := $(shell which gcc)
    CMAKE_CXX := $(shell which g++)
    LDFLAGS += -fopenmp
    LIBS = muparser.dll ./json/single_include/nlohmann/json.hpp
endif

all: $(LIBS)
	@$(CXX) $(CXXFLAGS) -O2 $(CPPFLAGS) -I./muparser/include -I./json/single_include/nlohmann $(SRCS) -o $(EXEC) -Wl,-rpath,. -L. -lmuparser

debug: $(LIBS)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I./muparser/include -I./json/single_include/nlohmann -DDEBUG $(SRCS) -o $(EXEC) -Wl,-rpath,. -L. -lmuparser

test: $(LIBS)
	@$(MAKE) --no-print-directory clean
	@cd ./test && mkdir output
	@mkdir plots
	@$(CXX) $(CXXFLAGS) -O2 $(CPPFLAGS) -I./muparser/include -I./json/single_include/nlohmann -I./test/ -DTEST $(SRCS_TEST) -o $(EXEC) -Wl,-rpath,. -L. -lmuparser
	@bash ./test/serial_run.sh
	@$(CXX) $(CXXFLAGS) -O2 $(CPPFLAGS) -I./muparser/include -I./json/single_include/nlohmann -I./test/ -DTEST -DPARALLEL $(SRCS_TEST) -o $(EXEC) -Wl,-rpath,. -L. -lmuparser
	@bash ./test/parallel_run.sh
	@$(RM) ./executable


libs:
	@echo "Cleaning previous libraries..."
	@$(MAKE) --no-print-directory allclean
	@echo "Cloning libraries' source files from github: json..."
	@git clone --depth 1 $(JSON_URL) $(JSON_DIR) > /dev/null 2>&1
	@echo "Cloning libraries' source files from github: muparser..."
	@git clone --depth 1 $(MUPARSER_URL) $(MUPARSER_DIR) > /dev/null 2>&1
	@cd $(MUPARSER_DIR) && \
	mkdir build > /dev/null 2>&1 && \
	cd build > /dev/null 2>&1 && \
	echo "Building shared muparser library..." && \
	CC=$(CMAKE_CC) CXX=$(CMAKE_CXX) cmake .. -DENABLE_SAMPLES=OFF -DENABLE_OPENMP=OFF > /dev/null 2>&1 && \
	$(MAKE) > /dev/null 2>&1
	@$(MAKE) --no-print-directory symlinks
	@echo "Libraries cloned and built successfully!"

symlinks:
ifeq ($(UNAME_S),Darwin)
	@echo "Creating symbolic link: libmuparser.dylib"
	@ln -sf ./muparser/build/libmuparser.dylib libmuparser.dylib
	@echo "Creating symbolic link: libmuparser.2.dylib"
	@ln -sf ./muparser/build/libmuparser.2.dylib libmuparser.2.dylib
else ifeq ($(UNAME_S),Linux)
	@echo "Creating symbolic link: libmuparser.so"
	@ln -sf ./muparser/build/libmuparser.so.2.3.5 libmuparser.so
	@echo "Creating symbolic link: libmuparser.so.2"
	@ln -sf ./muparser/build/libmuparser.so.2.3.5 libmuparser.so.2
else ifeq ($(UNAME_S),MINGW)
	@echo "Copying library: muparser.dll"
    @copy /Y .\muparser\build\muparser.dll muparser.dll
else ifeq ($(UNAME_S),Windows)
	@echo "Copying library: muparser.dll"
    @copy /Y .\muparser\build\muparser.dll muparser.dll
else
    $(error Unsupported operating system: $(UNAME_S))
endif
	
allclean:
	@cd ./test && $(RM) *.txt *.png 
	@$(RM) *.o *.a *.so *.so.2 *.dylib *.dll *.vtk *.txt *.png
	@$(RM) -f $(EXEC) 
	@$(RM) -rf $(DOCS)
	@$(RM) -rf $(JSON_DIR)
	@$(RM) -rf $(MUPARSER_DIR)

clean:
	@cd ./test && $(RM) -rf output
	@$(RM) -rf plots
	@$(RM) *.o *.a *.vtk *.txt *.png
	@$(RM) -f $(EXEC) 
	@$(RM) -rf $(DOCS)

doc:
	@$(DOXYGEN) $(DOXYFILE)

