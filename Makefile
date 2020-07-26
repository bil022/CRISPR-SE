CC=g++

FLAGS=-lpthread -std=gnu++11 -Wall
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	FLAGS+= -static
endif

SRC = pair.cc CREST.cc CrisprDB.cc
DEPS = CREST.h util.h
OBJ = pair.o CREST.o CrisprDB.o

se: $(OBJ)
	$(CC) -o $@ $^ $(FLAGS)

se-regex: $(OBJ)
	$(CC) -o $@ $^ $(FLAGS) -lpcre2-8

ref=../ref/hg38_chrM
test:
	./crest -p 4 -r $(ref)
	(cat $(ref).h; sort -T. -k3,3 -k4,4n $(ref).mm4) | samtools view -Sb - > $(ref).bam

all:
	@echo $(FLAGS)
