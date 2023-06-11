/* monolis_wrap_scalapack.h */
#ifndef MONOLIS_WRAP_SCALAPACK_H
#define MONOLIS_WRAP_SCALAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

void monolis_scalapack_gesvd_R(
  int      N_loc,
  int      M,
  double** A,
  double** S,
  double*  V,
  double** D,
  int      comm);

void monolis_scalapack_gesvd_R_c_main(
  int      N_loc,
  int      M,
  int      P,
  double*  A,
  double*  S,
  double*  V,
  double*  D,
  int      comm);

#ifdef __cplusplus
}
#endif

#endif
