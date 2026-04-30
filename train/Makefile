CC = gcc
CFLAGS = -O3 -march=native -fopenmp -ffast-math
LDFLAGS = -lm

TARGET = train_hmm
SRC = train_hmm.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean