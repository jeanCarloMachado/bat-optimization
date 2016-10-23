all: source

source:
	gcc src/main.c -lm -o bat
run: source
	./bat
clear:
	rm -rf bat
	rm -rf dump/*
