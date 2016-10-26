all: source

source:
	gcc src/main.c -lm -o bat
run: source
	./bat
	notify-send "Execution finished"
clear:
	rm -rf bat
	rm -rf dump/*
