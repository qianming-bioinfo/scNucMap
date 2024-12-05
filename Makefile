# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2  -pthread

# Directories
SRC_DIR = src
BIN_DIR = bin
PYTHON_DIR = python
R_DIR = r
SCRIPTS_DIR = scripts

# Source files
MAIN_SRCS = $(wildcard $(SRC_DIR)/calSDMat.cpp $(SRC_DIR)/countNuclFrag.cpp)
SUPPORT_SRCS = $(wildcard $(SRC_DIR)/myoperation.cpp)

# Object files
MAIN_OBJS = $(MAIN_SRCS:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)
SUPPORT_OBJS = $(SUPPORT_SRCS:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

# Executables
EXECUTABLES = $(BIN_DIR)/calSDMat $(BIN_DIR)/countNuclFrag

# Default target
all: $(EXECUTABLES) check-python-deps check-r-deps install-scripts-permission part-clean

# Rule for calSDMat executable
$(BIN_DIR)/calSDMat: $(BIN_DIR)/calSDMat.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# Rule for countNuclFrag executable
$(BIN_DIR)/countNuclFrag: $(BIN_DIR)/countNuclFrag.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compile main source files
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MF $(BIN_DIR)/$*.d -c $< -o $@

# Check and install Python dependencies
check-python-deps:
	@echo "Checking and installing Python dependencies..."
	pip3 install -r $(PYTHON_DIR)/requirements.txt

# Check and install R dependencies
check-r-deps:
	@echo "Checking and installing R dependencies..."
	Rscript $(R_DIR)/check_packages.R

# Install executable permissions for shell scripts
install-scripts-permission:
	@chmod +x $(SCRIPTS_DIR)/*

# Clean
part-clean:
	@rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.d

clean:
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.d $(EXECUTABLES)

# Include dependency files if they exist
-include $(BIN_DIR)/*.d

.PHONY: all check-python-deps check-r-deps clean install-scripts-permission
