#Makefile for roughly similarity search
#This file should be used with the file "topK_Rough.cpp"
CC=g++ -Wno-deprecated -DNDEBUG -O3
CFLAGS=-c

OFILES=time.o

all:rough_search

rebuild:clean all
	
clean:
	rm -f *.o

rough_search:${OFILES} topK_Rough.o
	${CC} -o $@ ${OFILES} topK_Rough.o

# ------------------------------------------
topK_Rough.o:topK_Rough.cpp
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}
	

# COMMON OBJECT FILES

time.o:Time.cpp Time.h
	${CC} ${CFLAGS} $< -o $@ ${INCLUDE}