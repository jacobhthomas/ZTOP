#include <cstdio>
#include <math.h>

#define nqvec 4

// Cuda Kernels
__global__ void fetch_idx
  (cudaTextureObject_t, int *, double *);
__global__ void texturevals
  (cudaTextureObject_t, double *);	
__global__ void polint4f
  (double *, double *, double *, double *);
__global__ void my_print
  (cudaTextureObject_t texObject, double * output) ;
__global__ void fetch_idx
  (cudaTextureObject_t texObject, int idx, double* out) ;
__global__ void texturevals
  (cudaTextureObject_t obj, double * output) ;

// extern functions called from fortran
extern "C" void polint4f_
  (double *, double *, double *, double *);
extern "C" void texturevals_kernel_
  (cudaTextureObject_t * obj, int size, double* out);
extern "C" void polint_wrapper_
  (double * , double *, double *, double *);
extern "C" void kernel_
  (cudaTextureObject_t *, int *);
extern "C" void print_kernel_
  (cudaTextureObject_t * obj, int * size) ;
extern "C" void fetch_
  (cudaTextureObject_t * obj, int idx, double* out) ;
extern "C" void texturevals_kernel_
  (cudaTextureObject_t * obj, int size, double* out) ;
extern "C" cudaTextureObject_t setup_
  (double * data, int * size);
extern "C" void free_cuda_memory_
  ();
extern "C" void partonx12_wrapper_
  (int *, int *, int *, int *, int *, int *, int *,
   double *, double *, double *,
   cudaTextureObject_t *, cudaTextureObject_t *, cudaTextureObject_t *,
   double *) ;

// Regular Cuda C definitions
void bin_search (double * data, double val, int * lo, int * hi);
void bin_search_check_x(int lo, int upper, int size, int * ret);
void bin_search_check_q(int lower, int upper, int size, int * ret) ;

// macros
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert
(cudaError_t code, const char *file, int line, bool abort=true) {
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

static __inline__ __device__ double
fetch_double(uint2 p){ return __hiloint2double(p.y, p.x);}
