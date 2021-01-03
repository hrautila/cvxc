
INC =  -Wall -I../../src/include -I$(HOME)/.local/include -O1 -g
LIB = -L$(HOME)/.local/lib

STATIC = -Wl,-Bstatic,$(LIB),-larmasd

all: floorplan.so

floorplan.o: floorplan.c
	$(CC) -c -fPIC $(INC) $(CFLAGS) -o $@ floorplan.c

floorplan.so: floorplan.o
	$(CC) -fPIC -shared $(LIB) -o $@ floorplan.o $(STATIC) -Wl,-Bdynamic -lm

