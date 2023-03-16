!> 非零構造決定モジュール（メイン関数群）
module mod_monolis_spmat_nonzero_pattern_util_wrap
  use mod_monolis_utils
  use mod_monolis_spmat_nonzero_pattern_util
  use iso_c_binding

  implicit none

contains

  subroutine monolis_get_CSC_format_c(NC, NR, NZ, index, item, indexR, itemR, permR) &
    & bind(c, name = "monolis_get_CSC_format")
    implicit none
    integer(c_int), intent(in), value :: NC, NR, NZ
    integer(c_int), intent(in), target :: index(0:NC)
    integer(c_int), intent(in), target :: item(NZ)
    integer(c_int), target :: indexR(0:NR)
    integer(c_int), target :: itemR(NZ)
    integer(c_int), target :: permR(NZ)
    integer(kint), pointer :: indexRt(:)
    integer(kint), pointer :: itemRt(:)
    integer(kint), pointer :: permRt(:)

    indexRt => indexR
    itemRt => itemR
    permRt => permR
    call monolis_get_CSC_format(NC, NR, NZ, index, item, indexRt, itemRt, permRt)
  end subroutine monolis_get_CSC_format_c

end module mod_monolis_spmat_nonzero_pattern_util_wrap
