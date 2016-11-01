#include <stdio.h>

int main()
{
    int *device;
    struct cudaDeviceProp *deviceProp;

    device = (int *) malloc(sizeof(int));
    deviceProp = (cudaDeviceProp *) malloc(sizeof(cudaDeviceProp));

    cudaGetDevice(device);
    cudaGetDeviceProperties(deviceProp, *device);


    printf("Name: %s\n", deviceProp->name);
    printf("Global memory: %u\n", deviceProp->totalGlobalMem);
}
