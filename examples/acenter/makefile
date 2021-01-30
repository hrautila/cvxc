
INC =  -Wall -I../../src/include -I$(HOME)/.local/include -O1 -g
LIB = -L$(HOME)/.local/lib
LIBCVXC := $(shell find ../../build -name libcvxc.a | head -1)
STATIC = -Wl,-Bstatic,$(LIB),-larmasd

all: acenter.so

acenter.o: acenter.c
	$(CC) -c -fPIC $(INC) $(CFLAGS) -o $@ acenter.c

acenter.so: acenter.o
	$(CC) -fPIC -shared -o $@ acenter.o $(LIBCVXC) $(STATIC) -Wl,-Bdynamic -lm

