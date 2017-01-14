CC:="gcc"
DEBUG="-g"
all: source
.PHONY: paper
.PHONY: tests

source:
	${CC} ${DEBUG} -c src/bat.c 
	${CC} ${DEBUG} -c src/common.c 
	${CC} ${DEBUG} bat.o common.o  src/main.c -lm -o bat 

run: source
	./bat
	./post-run

gpu:
	nvcc -c src/common.c -Wno-deprecated-gpu-targets 
	nvcc common.o src/main.cu -Wno-deprecated-gpu-targets -lm -o bat_gpu

run_gpu: gpu
	./bat_gpu
	./post-run

device_info:
	nvcc src/inspect_device.cu -o ls_device
	./ls_device

clear:
	rm -rf bat
	rm -rf dump/*
	rm -rf *.o

tests:
	${CC} ${DEBUG} -c src/common.c 
	${CC} ${DEBUG} -c src/bat.c 
	${CC} ${DEBUG} -c src/unity.c 
	${CC} ${DEBUG}  common.o bat.o unity.o src/tests.c -lm -o bat_tests 

run_tests: tests
	./bat_tests

paper:
	cd paper ; \
	pdflatex -interaction=nonstopmode paper.tex ; \
	mv paper.pdf  ../paper.pdf
