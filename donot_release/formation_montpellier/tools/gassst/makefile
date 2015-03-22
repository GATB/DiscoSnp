CC=g++
CFLAGS= -Wall -O3 -D_REENTRANT -msse2 -mssse3  -Wno-write-strings -DSSE=3 
LDFLAGS=-lm -Xlinker -lpthread -O3 -mssse3   
EXEC=Gassst gassst_to_sam
#SRC=$(wildcard *.cpp)
SRC=Alignment.cpp  code.cpp  Doublon.cpp  gapless.cpp Index.cpp  Pool.cpp  Stat.cpp Bank.cpp       display.cpp  filter.cpp   gassst.cpp   Hit.cpp misc.cpp   Seed.cpp  Server.cpp   withgap.cpp
OBJ= $(SRC:.cpp=.o)


#sans ssse3 pour amd
ifeq ($(amd),1)
	CFLAGS= -Wall -O3 -D_REENTRANT -msse2  -Wno-write-strings  -DSSE=2 
endif


ifeq ($(dbg),1)
	CFLAGS= -Wall -O0 -D_REENTRANT -msse2  -Wno-write-strings  -DSSE=2  -DDEBUG -g
	LDFLAGS+= -g
endif


ifeq ($(osx),1)
	CFLAGS+= -DOSX
endif
all: $(EXEC)

Gassst: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) 

gassst.o: gassst.cpp  code.h
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)


gassst_to_sam: gassst_to_sam.cpp
	$(CC)  $(CFLAGS) -o $@ $< 

.PHONY: clean

clean:
	rm -rf *.o