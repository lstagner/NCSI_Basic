static void HandleError(cudaError_t error, const char* file, int line){

  if( error != cudaSuccess){
    printf( "%s in %s at line %d\n",
	    cudaGetErrorString(error), file, line);
    exit( EXIT_FAILURE );
  }

} 

#define HANDLE_ERROR(error) (HandleError( error, __FILE__, __LINE__))
#define CHECK_KERNEL() (HandleError( cudaGetLastError(), __FILE__, __LINE__))
