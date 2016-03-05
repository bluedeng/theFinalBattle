#Makefile for similarity search in complex structure databases
CC=g++ -Wno-deprecated -DNDEBUG -O3
CFLAGS=-c

OFILES1 = time.o gram.o list.o filtb.o query.o db.o
OFILES2 = time.o gram.o list.o filtb.o query.o appro_db.o

all:topksearch approximateTopksearch roughsearch

rebuild:clean all
	
clean:
	rm -f *.o

topksearch:${OFILES1} topksearch.o
	${CC} -o $@ ${OFILES1} topksearch.o -lpthread

approximateTopksearch:${OFILES2} appro_topksearch.o
	${CC} -o $@ ${OFILES2} appro_topksearch.o -lpthread

roughsearch:time.o rough.o
	${CC} -o $@ time.o rough.o

# ------------------------------------------
topksearch.o:topkSearch.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}

appro_topksearch.o:Approximate_topkSearch.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	
rough.o:topK_Rough.cpp
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

appro_db.o:Approximate_SeqDB.cpp Approximate_SeqDB.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}