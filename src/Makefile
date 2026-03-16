CXX = g++
CXXFLAGS = -O3 -std=c++17 -march=native -ffast-math

LDFLAGS = -lm

TARGET = run_cern
SRC = run_cern.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
