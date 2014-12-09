CC=g++
KEYS=-std=c++11 -ldl -lpthread -o FASTA
SOURCE=FASTA.cpp

RELEASE: sqlite3.o
	$(CC) sqlite3.o $(KEYS) -O3 $(SOURCE) 

DEBUG: sqlite3.o
	$(CC) sqlite3.o $(KEYS) -ggdb $(SOURCE) 

sqlite3.o:
	gcc -c sqlite3.c -o sqlite3.o -lpthread -ldl

clean:
	rm -rf *.o FASTA