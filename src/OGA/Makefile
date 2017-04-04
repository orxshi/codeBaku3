CC = mpicxx
CFLAGS = -c -std=c++11 -fpic -O3
LDFLAGS = -shared
SOURCES=IBlank.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=libIBlank.so

all: $(TARGET)
	
$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean
clean:
	rm -rf $(OBJECTS) $(TARGET)
