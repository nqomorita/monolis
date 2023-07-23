/* monolis_spmat_copy.h */
#ifndef MONOLIS_SPMAT_COPY_H
#define MONOLIS_SPMAT_COPY_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 行列構造体のコピー（実数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列構造体のコピー（複素数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列の非零パターンのコピー（実数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列の非零パターンのコピー（複素数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列の非零パターンのコピー（実数型）
 * @param[in] NP 全計算点数
 * @param[in] NDOF 計算点が持つ自由度
 * @param[in] NZ 非零要素数
 * @param[in] NZU 非零要素数（上三角）
 * @param[in] NZL 非零要素数（下三角）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_val_R(
  int      NP,
  int      NDOF,
  int      NZ,
  int      NZU,
  int      NZL,
  MONOLIS_MAT_VAL_R* mat_in,
  MONOLIS_MAT_VAL_R* mat_out);

/**
 * @brief 行列の非零パターンのコピー（複素数型）
 * @param[in] NP 全計算点数
 * @param[in] NDOF 計算点が持つ自由度
 * @param[in] NZ 非零要素数
 * @param[in] NZU 非零要素数（上三角）
 * @param[in] NZL 非零要素数（下三角）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_val_C(
  int      NP,
  int      NDOF,
  int      NZ,
  int      NZU,
  int      NZL,
  MONOLIS_MAT_VAL_C* mat_in,
  MONOLIS_MAT_VAL_C* mat_out);

/**
 * @brief 行列の非零パターンのコピー（CSR）
 * @param[in] NP 全計算点数
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_CSR(
  int      NP,
  MONOLIS_MAT_CSR* mat_in,
  MONOLIS_MAT_CSR* mat_out);

/**
 * @brief 行列の非零パターンのコピー（CSC）
 * @param[in] NP 全計算点数
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_CSC(
  int      NP,
  MONOLIS_MAT_CSC* mat_in,
  MONOLIS_MAT_CSC* mat_out);

/**
 * @brief 行列の非零パターンのコピー（SCSR）
 * @param[in] NP 全計算点数
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_nonzero_pattern_SCSR(
  int      NP,
  MONOLIS_MAT_SEPARATED_CSR* mat_in,
  MONOLIS_MAT_SEPARATED_CSR* mat_out);

/**
 * @brief 行列値のコピー（実数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（行列、実数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_matrix_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（解ベクトル、実数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_solution_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（右辺ベクトル、実数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_rhs_R(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（複素数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（行列、複素数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_matrix_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（解ベクトル、複素数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_solution_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のコピー（右辺ベクトル、複素数型）
 * @param[in] mat_in monolis 構造体（入力)
 * @param[inout] mat_out monolis 構造体（出力）
 * @ingroup matrix_copy
 */
void monolis_copy_mat_value_rhs_C(
  MONOLIS* mat_in,
  MONOLIS* mat_out);

/**
 * @brief 行列値のゼロ初期化（実数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_R(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（行列、実数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_matrix_R(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（解ベクトル、実数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_solution_R(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（右辺ベクトル、実数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_rhs_R(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（複素数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_C(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（行列、複素数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_matrix_C(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（解ベクトル、複素数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_solution_C(
  MONOLIS* mat);

/**
 * @brief 行列値のゼロ初期化（右辺ベクトル、複素数型）
 * @param[inout] mat monolis 構造体（入力)
 * @ingroup matrix_copy
 */
void monolis_clear_mat_value_rhs_C(
  MONOLIS* mat);

#ifdef __cplusplus
}
#endif

#endif
