CC=gcc
OPTS=-c -fPIC
LOADER=gcc

SRC=leven.c hashmap.c

leven.so: $(SRC)
	R CMD SHLIB $(SRC)

clean:
	rm *.o *.so
