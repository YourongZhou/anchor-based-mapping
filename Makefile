CXX=clang++
CXXFLAGS=-std=c++11 -fPIC -O2 -Wall
LDFLAGS=
SEQAN_INC=`pkg-config --cflags --libs seqan`  # 如果没有 pkg-config，可手动设置 -I/path/to/seqan/include

INCLUDE_DIR=include
SRC_DIR=src
BUILD_DIR=build

LIBS=$(BUILD_DIR)/libpart1.so $(BUILD_DIR)/libpart2.so $(BUILD_DIR)/libpart3.so

all: dirs $(BUILD_DIR)/main_test

dirs:
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/libpart1.so: $(SRC_DIR)/part1.cpp $(INCLUDE_DIR)/part1.h $(INCLUDE_DIR)/anchor_common.h
	$(CXX) $(CXXFLAGS) -shared -o $@ $(SRC_DIR)/part1.cpp $(SEQAN_INC)

$(BUILD_DIR)/libpart2.so: $(SRC_DIR)/part2.cpp $(INCLUDE_DIR)/part2.h $(INCLUDE_DIR)/anchor_common.h
	$(CXX) $(CXXFLAGS) -shared -o $@ $(SRC_DIR)/part2.cpp $(SEQAN_INC)

$(BUILD_DIR)/libpart3.so: $(SRC_DIR)/part3.cpp $(INCLUDE_DIR)/part3.h $(INCLUDE_DIR)/anchor_common.h
	$(CXX) $(CXXFLAGS) -shared -o $@ $(SRC_DIR)/part3.cpp $(SEQAN_INC)

$(BUILD_DIR)/main_test: $(SRC_DIR)/main_test.cpp $(BUILD_DIR)/libpart1.so $(BUILD_DIR)/libpart2.so $(BUILD_DIR)/libpart3.so
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/main_test.cpp -L$(BUILD_DIR) -lpart1 -lpart2 -lpart3 $(SEQAN_INC) -Wl,-rpath,'$$ORIGIN'

clean:
	rm -rf $(BUILD_DIR) anchors.fasta *.tmp.txt

.PHONY: all clean dirs
