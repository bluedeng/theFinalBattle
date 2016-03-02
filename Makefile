#Makefile for similarity search in complex structure databases
CC=g++ -Wno-deprecated -DNDEBUG -O3
CFLAGS=-c

OFILES=time.o gram.o list.o filtb.o query.o db.o

all:psearch

rebuild:clean all
	
clean:
	rm -f *.o

psearch:${OFILES} topksearch.o
	${CC} -o $@ ${OFILES} psearch.o -lpthread

# ------------------------------------------
topksearch.o:topkSearch.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	

# COMMON OBJECT FILES

time.o:Time.cpp Time.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
gram.o:Gram.cpp Gram.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

list.o:GramList.cpp GramList.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
filtb.o:CountFilter.cpp CountFilter.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

query.o:Query.cpp Query.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

db.o:SeqDB.cpp SeqDB.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
