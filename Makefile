CC:="gcc"
DEBUG="-g"
all: cpu gpu
.PHONY: paper
.PHONY: tests

test:
	${CC} ${DEBUG} -c src/bat/mersenne.c
	${CC} ${DEBUG} -c src/bat/common.c
	${CC} ${DEBUG} -c src/bat/cpu.c
	${CC} ${DEBUG} -c src/unity.c
	${CC} ${DEBUG} cpu.o mersenne.o common.o unity.o src/test_cpu.c -lm -o bat_tests
	./bat_tests

cpu:
	${CC} ${DEBUG} -c src/bat/mersenne.c
	${CC} ${DEBUG} -c src/bat/common.c
	${CC} ${DEBUG} -c src/bat/cpu.c
	${CC} ${DEBUG} cpu.o common.o mersenne.o src/main.c -lm -o bat

run: cpu
	./bat

gpu:
	nvcc -c src/common.c -Wno-deprecated-gpu-targets
	nvcc common.o src/main.cu -Wno-deprecated-gpu-targets -lm -o bat_gpu

run_gpu: gpu
	./bat_gpu

device_info:
	nvcc src/inspect_device.cu -o ls_device
	./ls_device

clear:
	rm -rf bat
	rm -rf dump/*
	rm -rf *.o



paper:
	cd paper ; \
	pdflatex paper.tex ; \
	mv paper.pdf  ../paper.pdf
