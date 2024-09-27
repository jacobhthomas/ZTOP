#include<stdio.h>

__global__ void parentkern(double v[8])
{
  v[0] = 3.5;
}


extern "C" void simple_
(int * i1, double * v1, double vout[8])
{
  vout[2] = 2.3;
  double* device_v;
  cudaMalloc((void**)&device_v, sizeof(double) * 8);
  
  parentkern<<<1,1>>>(device_v);
  cudaDeviceSynchronize();
  cudaMemcpy(vout, device_v, sizeof(double)*8, cudaMemcpyDeviceToHost);

  cudaFree(device_v);

  //  printf("%d\n",cudaGetLastError());
  printf("%s",cudaGetErrorString(cudaGetLastError()));
}

