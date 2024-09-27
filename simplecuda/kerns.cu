#include<stdio.h>

__device__ int ii1;

__global__ void childkern(int i1, double v1, double vout[8])
{
   int i = threadIdx.x;
   vout[i] = i;
}


__global__ void parentkern(int* i1, double* v1, double vout[8])
{
   ii1 = *i1;
   double vv1 = *v1;

   ii1 += 3;
//   childkern<<<1,8>>>(ii1, vv1, vout);
   vout[0] = 2.3;
      
}


extern "C" void kerns_
(int * i1, double * v1, double vout[8])
{
//  printf("%f\n",*vout[0]);
//  printf("%f\n",*vout[1]);

  int * device_i1;
  double * device_v1; // device vars
  double * device_vout ;

  // allocate memory on device (GPU)
  cudaMalloc((void**)&device_i1 , sizeof(int));
  cudaMalloc((void**)&device_v1 , sizeof(double));
  cudaMalloc((void**)&device_vout, sizeof(double) * 8);

  // copy sent data to GPU
  cudaMemcpy(device_i1 , i1, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(device_v1, v1, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_vout, vout, sizeof(double)*8, cudaMemcpyHostToDevice);

  // call on GPU
  parentkern<<<1,1>>>(device_i1, device_v1, device_vout);

  // copy value of Y to cpu
  cudaMemcpy(vout, device_vout, sizeof(double)*8, cudaMemcpyDeviceToHost);
  // free memory
  cudaFree(device_vout);
  cudaFree(device_v1);
  cudaFree(device_i1);
}
