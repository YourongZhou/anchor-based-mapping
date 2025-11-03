# Compiler
CXX      = clang++
CXXFLAGS = -std=c++17 -O3 -fPIC -Iinclude
# CXXFLAGS = -std=c++17 -g -O0 -fPIC -Iinclude // debug

# Source and object files
SRC = src/anchor_gen.cpp \
      src/candidate_retriever.cpp \
      src/index_query.cpp \
      src/metrics.cpp \
      src/rng.cpp \
      src/tools.cpp \
      main.cpp
#       src/test.cpp \

OBJ = $(SRC:.cpp=.o)

# Output binary
TARGET = anchor_mapping

# Default rule
all: $(TARGET)

# Link final executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

# Compile each .cpp into .o
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f src/*.o $(TARGET)
