CC = gcc
CFLAGS = -Wall -O3 -shared -fPIC -lm

C_SRC = c_src
LIB_DIR = lib

TARGETS = $(LIB_DIR)/rawhash_events.so $(LIB_DIR)/rawhash2_events.so

all: $(LIB_DIR) $(TARGETS)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

$(LIB_DIR)/rawhash_events.so: $(C_SRC)/rawhash_events.c | $(LIB_DIR)
	$(CC) $(CFLAGS) -o $@ $<

$(LIB_DIR)/rawhash2_events.so: $(C_SRC)/rawhash2_events.c | $(LIB_DIR)
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -rf $(LIB_DIR)

.PHONY: all clean
