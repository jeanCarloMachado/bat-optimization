all: source

source:
	gcc -c src/mersenne.c
	gcc mersenne.o src/main.c  -lm -o bat 
run: source
	./bat
	./post-run


source_gpu:
	nvcc -c src/mersenne.c
	nvcc mersenne.o src/main.cu -Wno-deprecated-gpu-targets -lm -o bat_gpu

run_gpu: source_gpu
	./bat_gpu

device_info:
	nvcc src/inspect_device.cu -o ls_device
	./ls_device

clear:
	rm -rf bat
	rm -rf dump/*
