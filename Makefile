# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11  -O3 -march=native -fopenmp -Iinc -I/usr/include/eigen3 -g

# Directories
SRC_DIR = src
INC_DIR = inc
OBJ_DIR = bin
MN_DIR = main

# Target executable
#TARGET = single_e-

# Default program to run (can be overridden: make run PROGRAM=...)
PROGRAM ?= $(TARGET)

# Source files and object files
SRCS = $(SRC_DIR)/functions.cpp $(SRC_DIR)/iofile.cpp main/$(TARGET).cpp
OBJS = $(SRCS:.cpp=.o)

# Default target: build and run
all: run

# Build only
build: $(TARGET)

# Run the program (after building)
run: build
	@if [ -x $(PROGRAM) ]; then \
		echo "Running $(PROGRAM)..."; \
		./$(PROGRAM); \
	else \
		echo "Error: Program '$(PROGRAM)' not found or not executable."; \
		exit 1; \
	fi

# Link object files into final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile .cpp to .o
%.o: %.cpp $(INC_DIR)/functions.h $(INC_DIR)/iofile.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(SRC_DIR)/*.o *.o $(TARGET) $(PROGRAM) $(MN_DIR)/*.o

# Declare phony targets
.PHONY: all build run clean
