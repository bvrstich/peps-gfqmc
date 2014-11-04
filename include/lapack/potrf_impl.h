#ifndef __BTAS_LAPACK_POTRF_IMPL_H
#define __BTAS_LAPACK_POTRF_IMPL_H 1

#include <lapack/types.h>

namespace btas {
   namespace lapack {

      template<typename T>
         void potrf (
               const int& order,
               const char &uplo,
               const size_t& N,
               T* A,
               const size_t& ldA)
         {
            BTAS_LAPACK_ASSERT(false, "potrf must be specialized.");
         }

      inline void potrf (
            const int& order,
            const char &uplo,
            const size_t& N,
            float* A,
            const size_t& ldA)
      {
         LAPACKE_spotrf(order,uplo, N, A, ldA);
      }

      inline void potrf (
            const int& order,
            const char &uplo,
            const size_t& N,
            double* A,
            const size_t& ldA)
      {
         LAPACKE_dpotrf(order,uplo, N, A, ldA);

      }

      inline void potrf (
            const int& order,
            const char &uplo,
            const size_t& N,
            std::complex<float>* A,
            const size_t& ldA)
      {
         LAPACKE_cpotrf(order,uplo, N, A, ldA);
      }

      inline void potrf (
            const int& order,
            const char &uplo,
            const size_t& N,
            std::complex<double>* A,
            const size_t& ldA) {
         LAPACKE_zpotrf(order, uplo,N, A, ldA);
      }

   } // namespace lapack
} // namespace btas

#endif // __BTAS_LAPACK_POTRF_IMPL_H
