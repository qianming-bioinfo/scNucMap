# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -pthread

# Directories
SRC_DIR = src
BIN_DIR = bin
PYTHON_DIR = python
R_DIR = r
SCRIPTS_DIR = scripts

# Source files
# MAIN_SRCS = $(wildcard $(SRC_DIR)/cal_SD_and_frag_count_lite.cpp \
# 					   $(SRC_DIR)/cal_SD_and_frag_count_standard.cpp \
# 					   $(SRC_DIR)/calSDMat_lite.cpp \
# 					   $(SRC_DIR)/calSDMat_standard.cpp \
# 					   $(SRC_DIR)/countNuclFrag_lite.cpp \
# 					   $(SRC_DIR)/countNuclFrag_standard.cpp)

MAIN_SRCS = $(wildcard $(SRC_DIR)/cal_SD_and_frag_count_lite.cpp \
					   $(SRC_DIR)/cal_SD_and_frag_count_standard.cpp \
					   $(SRC_DIR)/calSDMat_lite.cpp \
					   $(SRC_DIR)/calSDMat_standard.cpp \
					   $(SRC_DIR)/countNuclFrag_standard.cpp)

SUPPORT_SRCS = $(wildcard $(SRC_DIR)/myoperation.cpp)

# Object files
MAIN_OBJS = $(MAIN_SRCS:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)
SUPPORT_OBJS = $(SUPPORT_SRCS:$(SRC_DIR)/%.cpp=$(BIN_DIR)/%.o)

# Executables
EXECUTABLES = \
	$(BIN_DIR)/cal_SD_and_frag_count_lite \
	$(BIN_DIR)/cal_SD_and_frag_count_standard \
	$(BIN_DIR)/calSDMat_lite \
	$(BIN_DIR)/calSDMat_standard \
	$(BIN_DIR)/countNuclFrag_standard \
	# $(BIN_DIR)/countNuclFrag_lite


# Default target
all: $(EXECUTABLES) check-python-deps check-r-deps install-scripts-permission part-clean

# Rules for executables
$(BIN_DIR)/cal_SD_and_frag_count_lite: $(BIN_DIR)/cal_SD_and_frag_count_lite.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/cal_SD_and_frag_count_standard: $(BIN_DIR)/cal_SD_and_frag_count_standard.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/calSDMat_lite: $(BIN_DIR)/calSDMat_lite.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/calSDMat_standard: $(BIN_DIR)/calSDMat_standard.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# $(BIN_DIR)/countNuclFrag_lite: $(BIN_DIR)/countNuclFrag_lite.o $(BIN_DIR)/myoperation.o
# 	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/countNuclFrag_standard: $(BIN_DIR)/countNuclFrag_standard.o $(BIN_DIR)/myoperation.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compile source files
$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MF $(BIN_DIR)/$*.d -c $< -o $@

# Python dependencies
check-python-deps:
	@echo "Checking and installing Python dependencies..."
	pip3 install -r $(PYTHON_DIR)/requirements.txt

# R dependencies
check-r-deps:
	@echo "Checking and installing R dependencies..."
	Rscript $(R_DIR)/check_packages.R

# Set shell script executable permissions
install-scripts-permission:
	@chmod +x $(SCRIPTS_DIR)/*

# Clean only object and dependency files
part-clean:
	@rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.d

# Full clean
clean:
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/*.d $(EXECUTABLES)

# Install all (optional alias)
install: all

# Include dependency files
-include $(BIN_DIR)/*.d

.PHONY: all clean part-clean check-python-deps check-r-deps install-scripts-permission install
