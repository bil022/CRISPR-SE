FLAGS=-lpthread -Wall
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	FLAGS+= -static -lrt
endif

all:
	g++ *.cc -o se $(FLAGS)
debug:
	g++ *.cc -o ../seD -DDEBUG -g -pg $(FLAGS)
index:
	./se --index -r ecoli
build:
	./se --build -r ecoli
gz:
	echo tar -cvzf CRISPR-SE.tgz CRISPR-SE
	echo "pwd | mail -s CRISPR-SE -a CRISPR-SE.tgz bil022@ucsd.edu"
