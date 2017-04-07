CC:="gcc"
#perfromance flags
#FLAGS="-O3"
#debug flags
FLAGS="-g"
all: cpu gpu
.PHONY: paper
.PHONY: tests

test:
	${CC} ${FLAGS} -c src/bat/mersenne.c
	${CC} ${FLAGS} -c src/bat/common.c
	${CC} ${FLAGS} -c src/bat/cpu.c
	${CC} ${FLAGS} -c src/unity.c
	${CC} ${FLAGS} cpu.o mersenne.o common.o unity.o src/test_cpu.c -lm -o bat_tests
	./bat_tests

cpu: clear
	${CC} ${FLAGS} -c src/bat/mersenne.c
	${CC} ${FLAGS} -c src/bat/common.c
	${CC} ${FLAGS} -c src/bat/cpu.c
	${CC} ${FLAGS} cpu.o common.o mersenne.o src/main.c -lm -o bat

run: cpu
	./bat

gpu: clear
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
