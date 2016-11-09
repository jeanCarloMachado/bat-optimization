all: source
.PHONY: paper

source:
	gcc -c src/mersenne.c
	gcc -c src/common.c 
	gcc mersenne.o common.o src/main.c -lm -o bat 

run: source
	./bat
	./post-run


gpu:
	nvcc -c src/mersenne.c -Wno-deprecated-gpu-targets 
	nvcc -c src/common.c -Wno-deprecated-gpu-targets 
	nvcc common.o mersenne.o src/main.cu -Wno-deprecated-gpu-targets -lm -o bat_gpu

run_gpu: gpu
	./bat_gpu
	./post-run

device_info:
	nvcc src/inspect_device.cu -o ls_device
	./ls_device

clear:
	rm -rf bat
	rm -rf dump/*

paper:
	cd paper ; \
	pdflatex paper.tex
	mv paper/paper.pdf  paper.pdf
