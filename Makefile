all: source

source:
	gcc src/main.c -lm -o bat
clear:
	rm -rf bat
