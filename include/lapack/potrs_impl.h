#ifndef __BTAS_LAPACK_POTRS_IMPL_H
#define __BTAS_LAPACK_POTRS_IMPL_H 1

#include <lapack/types.h>

namespace btas {
   namespace lapack {

      template<typename T>
         void potrs (
               const int& order,
               const char &uplo,
               const size_t& N,
               const size_t& nrs,
               const T* A,
               const size_t& ldA,
               T* b,
               const size_t& ldb)
         {
            BTAS_LAPACK_ASSERT(false, "potrs must be specialized.");
         }

      inline void potrs (
            const int& order,
            const char &uplo,
            const size_t& N,
            const size_t& nrs,
            const float* A,
            const size_t& ldA,
            float* b,
            const size_t& ldb)
      {
         LAPACKE_spotrs(order,uplo, N,nrs, A, ldA,b,ldb);
      }

      inline void potrs (
            const int& order,
            const char &uplo,
            const size_t& N,
            const size_t& nrs,
            const double* A,
            const size_t& ldA,
            double* b,
            const size_t& ldb)
      {
         LAPACKE_dpotrs(order,uplo, N,nrs, A, ldA,b,ldb);

      }

      inline void potrs (
            const int& order,
            const char &uplo,
            const size_t& N,
            const size_t& nrs,
            const std::complex<float>* A,
            const size_t& ldA,
            std::complex<float>* b,
            const size_t& ldb)
      {
         LAPACKE_cpotrs(order,uplo, N, nrs, A, ldA, b, ldb);
      }

      inline void potrs (
            const int& order,
            const char &uplo,
            const size_t& N,
            const size_t& nrs,
            const std::complex<double>* A,
            const size_t& ldA,
            std::complex<double>* b,
            const size_t& ldb)
      {
         LAPACKE_zpotrs(order, uplo,N,nrs, A, ldA,b,ldb);
      }

   } // namespace lapack
} // namespace btas

#endif // __BTAS_LAPACK_POTRS_IMPL_H
