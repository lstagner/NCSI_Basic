// Macro for compiling on host and device
#ifdef __CUDACC__
#define CUDA_MEMBER __host__ __device__
#define HOST_MEMBER __host__
#define DEVICE_MEMBER __device__
#else
#define CUDA_MEMBER
#define HOST_MEMBER
#define DEVICE_MEMBER
#endif // __CUDACC__

