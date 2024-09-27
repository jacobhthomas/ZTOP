#include "defns.cuh"

__global__ void polint4f (double XA[4], double YA[4], double * X, double* Y) {
  double H1,H2,H3,H4,W,DEN,D1,C1,D2,C2,D3,C3,CD1,CC1,CD2,CC2,DD1,DC1;
  H1 = XA[0] - *X;
  H2 = XA[1] - *X;
  H3 = XA[2] - *X;
  H4 = XA[3] - *X;

  W   = YA[1] - YA[0];
  DEN = W / (H1 - H2);
  D1  = H2 * DEN;
  C1  = H1 * DEN;

  W   = YA[2] - YA[1];
  DEN = W/(H2-H3);
  D2  = H3*DEN;
  C2  = H2*DEN;

  W   = YA[3]-YA[2];
  DEN = W/(H3-H4);
  D3  = H4*DEN;
  C3  = H3*DEN;

  W   = C2-D1;
  DEN = W/(H1-H3);
  CD1 = H3*DEN;
  CC1 = H1*DEN;

  W   = C3-D2;
  DEN = W/(H2-H4);
  CD2 = H4*DEN;
  CC2 = H2*DEN;

  W   = CC2-CD1;
  DEN = W/(H1-H4);
  DD1 = H4*DEN;
  DC1 = H1*DEN;

  if(H3 + H4 < 0.0)
    *Y = YA[3] + D3 + CD2 + DD1;
  else if(H2 + H3 < 0.0)
    *Y = YA[2] + D2 + CD1 + DC1;
  else if(H1+H2 < 0.0)
    *Y = YA[1] + C2 + CD1 + DC1;
  else
    *Y = YA[0] + C1 + CC1 + DC1;
}

// Global device variables
double *device_XA, *device_YA, *device_X, *device_Y;

// gpu allocation
extern "C" void allocate_memory_(double XA[4], double YA[4], double* XV, double* YV) {
  cudaMalloc((void**)device_X , sizeof(double));
  cudaMalloc((void**)device_Y , sizeof(double));
  cudaMalloc((void**)device_XA, sizeof(double) * 4);
  cudaMalloc((void**)device_YA, sizeof(double) * 4);
}

// copy sent data to GPU
extern "C" void copy_memory_(double XA[4], double YA[4], double* XV, double* YV) {
  cudaMemcpy(device_X , XV, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_XA, XA, sizeof(double) * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(device_YA, YA, sizeof(double) * 4, cudaMemcpyHostToDevice);
}

// gpu deallocation
extern "C" void free_memory_(double XA[4], double YA[4], double* XV, double* YV) {
  cudaFree(device_X);
  cudaFree(device_Y);
  cudaFree(device_XA);
  cudaFree(device_YA);
}

extern "C" void polint_wrapper_(double XA[4], double YA[4], double* XV, double* YV) {
  double * device_XA, * device_YA; // device vars
  double * device_X , * device_Y ;

  // allocate memory on device (GPU)
  cudaMalloc((void**)&device_Y , sizeof(double));
  cudaMalloc((void**)&device_X , sizeof(double));
  cudaMalloc((void**)&device_XA, sizeof(double) * 4);
  cudaMalloc((void**)&device_YA, sizeof(double) * 4);

  // copy sent data to GPU
  cudaMemcpy(device_X , XV, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(device_XA, XA, sizeof(double) * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(device_YA, YA, sizeof(double) * 4, cudaMemcpyHostToDevice);

  // call on GPU
  polint4f<<<1,1>>>(device_XA, device_YA, device_X, device_Y);

  // copy value of Y to cpu
  cudaMemcpy(YV, device_Y, sizeof(double), cudaMemcpyDeviceToHost);

  // free memory
  cudaFree(device_X);
  cudaFree(device_Y);
  cudaFree(device_XA);
  cudaFree(device_YA);

}

extern "C" void polint_wrapper_modified_(double XA[4], double YA[4], double* XV, double* YV) {

  polint4f<<<1,1>>>(device_XA, device_YA, device_X, device_Y);

  cudaMemcpy(YV, device_Y, sizeof(double), cudaMemcpyDeviceToHost);

}