#include "defns.cuh"

__global__ void my_print(cudaTextureObject_t texObject, double * output) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  uint2 rval = tex1Dfetch<uint2>(texObject, idx);
  double dval = fetch_double(rval);
  __syncthreads();
  printf("TEX: %llu IDX: %d VAL: %0.16f\n", texObject, idx, dval);
  output[idx] = dval;
}

extern "C" void print_kernel_(cudaTextureObject_t * obj, int * size) {
  // Kernel for calling my_print from Fortran
  double * d_out, * out;
  cudaMalloc((void**)&d_out, sizeof(double) * (*size));
  out = (double*)malloc(sizeof(double) * (*size));
  
  // calls my_print from fortran
  my_print<<<*size,1>>>(*obj,d_out);
  cudaMemcpy(out,d_out, sizeof(double) * (*size), cudaMemcpyDeviceToHost);
  cudaFree(d_out);
}

__global__ void fetch_idx(cudaTextureObject_t texObject, int idx, double* out) {
  // Same as my_print, but yields a single value in device variable out
  uint2 rval = tex1Dfetch<uint2>(texObject, idx);
  *out = fetch_double(rval);
}

extern "C" void fetch_ (cudaTextureObject_t * obj, int idx, double* out) {
  /* given a texture object obj, an index idx, put
   obj[idx] in the single double value out */

  double * d_out;
  cudaMalloc((void**)&d_out, sizeof(double));
  fetch_idx<<<1,1>>>(*obj,idx,d_out);
  cudaMemcpy(out,d_out, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_out);
}

__global__ void texturevals (cudaTextureObject_t obj, double * output) {
  // fetches values in texture obj into gpu array output
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  uint2 rval = tex1Dfetch<uint2>(obj, idx);
  output[idx] = fetch_double(rval);
  __syncthreads();
 }

// kernel for texturevals
extern "C" void texturevals_kernel_
(cudaTextureObject_t * obj, int size, double* out) {
  
  double * d_out;
  cudaMalloc((void**)&d_out, size * sizeof(double));

  texturevals<<<size,1>>>(*obj,d_out);
  cudaMemcpy(out,d_out, size * sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(d_out);
}

extern "C" cudaTextureObject_t setup_(double * data, int * size){
  // initialize device data
  
  double* d_data;
  cudaMalloc((void**)&d_data,(*size)*sizeof(double));
  cudaMemcpy(d_data, data, (*size)*sizeof(double), cudaMemcpyHostToDevice);

  cudaTextureDesc td;
  memset(&td, 0, sizeof(td));
  td.normalizedCoords = 0;
  td.addressMode[0] = cudaAddressModeBorder;
  td.readMode = cudaReadModeElementType;

  struct cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeLinear;
  resDesc.res.linear.devPtr = d_data;
  resDesc.res.linear.sizeInBytes = *size*sizeof(double);
  resDesc.res.linear.desc.f = cudaChannelFormatKindUnsigned;
  resDesc.res.linear.desc.x = 32;
  resDesc.res.linear.desc.y = 32;

  cudaTextureObject_t texObject;
  gpuErrchk(cudaCreateTextureObject(&texObject, &resDesc, &td, NULL));

  // printf("LEAVING SETUP\n");  
  // printf("TEXID: %llu\n", texObject);
  // free(data);
  return texObject;
}

extern "C" void free_cuda_memory_() {
  cudaDeviceReset();
}