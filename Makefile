all: source

source:
	gcc src/main.c -lm -o bat
run: source
	./bat | cut -d " " -f3 | get-line 1 | mycopy
clear:
	rm -rf bat
	rm -rf dump/*
