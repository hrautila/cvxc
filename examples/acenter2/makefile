
INC =  -Wall -I../../src/include -I$(HOME)/.local/include -O1 -g
LIB = -L$(HOME)/.local/lib
LIBCVXC := $(shell find ../../build -name libcvxc.a | head -1)
STATIC = -Wl,-Bstatic,$(LIB),-larmasd

all: acenter_eq.so acenter_data

acenter_eq.o: acenter_eq.c
	$(CC) -c -fPIC $(INC) $(CFLAGS) -o $@ acenter_eq.c

acenter_eq.so: acenter_eq.o
	$(CC) -fPIC -shared $(LIB) -o $@ acenter_eq.o $(LIBCVXC) $(STATIC) -Wl,-Bdynamic -lm

acenter_data: acenter_data.c
	$(CC) -o $@ $(INC) acenter_data.c $(LIBCVXC) $(STATIC) -Wl,-Bdynamic -lm
