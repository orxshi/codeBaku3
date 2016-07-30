OBJDIR = obj
INCDIR = inc
CC = g++
CFLAGS = -c -std=c++11 -fpic -O3 -MMD -I$(INCDIR) -larmadillo
LDFLAGS = -shared
SRC := $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,obj/%.o,$(SRC)) 
TARGET=libcb.so

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $@

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@
  
-include $(ODBJIR)/*.d
	
.PHONY: clean
clean:
	rm -rf $(OBJ) $(TARGET)
