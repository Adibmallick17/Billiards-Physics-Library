CC = clang
CFLAGS = -Wall -std=c99 -pedantic
LDFLAGS = -L. -lm -lphylib

all: phylib

clean:
	rm -f *.o *.so phylib

libphylib.so: phylib.o
	$(CC) -shared -o libphylib.so phylib.o

phylib.o: phylib.c phylib.h
	$(CC) $(CFLAGS) -fPIC -c phylib.c -o phylib.o

A1test1.o: A1test1.c phylib.h
	$(CC) $(CFLAGS) -c A1test1.c -o A1test1.o

phylib: A1test1.o libphylib.so
	$(CC) A1test1.o $(LDFLAGS) -Wl,-rpath,. -o phylib
	




