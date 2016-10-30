all: source

source:
	gcc mersenne.o src/main.c  -lm -o bat 
run: source
	./bat
	./post-run
clear:
	rm -rf bat
	rm -rf dump/*
