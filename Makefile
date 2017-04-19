#perfromance flags
#FLAGS="-O3"
#debug flags
CC:="gcc"
GC:="nvcc"
FLAGS="-g"
GPU_FLAGS="-Wno-deprecated-gpu-targets"
all: cpu gpu
.PHONY: paper
.PHONY: tests

test:
	${CC} ${FLAGS} -c src/bat/cpu.c
	${CC} ${FLAGS} -c src/unity.c
	${CC} ${FLAGS} cpu.o unity.o src/test_cpu.c -lm -o bat_tests
	./bat_tests

cpu: clear
	${CC} ${FLAGS} -c src/bat/cpu.c
	${CC} ${FLAGS} cpu.o src/main.c -lm -o bat

run: cpu
	./bat

gpu: clear
	${GC} ${GPU_FLAGS} -x cu -c src/bat/gpu.cu
	${GC} ${GPU_FLAGS} -x cu -c src/main.c
	${GC} ${GPU_FLAGS} gpu.o main.o -lm -o bat_gpu
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
