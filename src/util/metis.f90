module mod_monolis_metis
  use, intrinsic :: iso_c_binding

  implicit none

  integer, parameter :: idx_t  = c_int
  integer, parameter :: real_t = c_float

  interface
    function metis_nodeND_f(nvtxs, xadj, adjncy, vwgt, options, perm, iperm) &
      result(ierr) bind(C, name="METIS_NodeND")
      import
      integer(idx_t)             :: ierr
      integer(idx_t), intent(in) :: nvtxs
      type(c_ptr), value         :: xadj,adjncy
      type(c_ptr), value         :: vwgt
      type(c_ptr), value         :: options
      type(c_ptr), value         :: perm
      type(c_ptr), value         :: iperm
    end function metis_nodeND_f
  end interface

contains

  function monolis_metis_nodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm) result(ierr)
    integer(idx_t) :: ierr
    integer(idx_t) :: nvtxs
    integer(idx_t), pointer :: xadj(:)
    integer(idx_t), pointer :: adjncy(:)
    integer(idx_t), pointer :: options(:)
    integer(idx_t), pointer :: vwgt(:)
    integer(idx_t), pointer :: perm(:)
    integer(idx_t), pointer :: iperm(:)

    type(c_ptr) :: xadj_ptr, adjncy_ptr, vwgt_ptr
    type(c_ptr) :: options_ptr
    type(c_ptr) :: perm_ptr
    type(c_ptr) :: iperm_ptr

    xadj_ptr    = c_null_ptr
    adjncy_ptr  = c_null_ptr
    vwgt_ptr    = c_null_ptr
    options_ptr = c_null_ptr
    perm_ptr    = c_null_ptr
    iperm_ptr   = c_null_ptr

    if(associated(xadj   )) xadj_ptr    = c_loc(xadj(1))
    if(associated(adjncy )) adjncy_ptr  = c_loc(adjncy(1))
    if(associated(vwgt   )) vwgt_ptr    = c_loc(vwgt(1))
    if(associated(options)) options_ptr = c_loc(options(1))
    if(associated(perm   )) perm_ptr    = c_loc(perm(1))
    if(associated(iperm  )) iperm_ptr   = c_loc(iperm(1))

#ifdef WITHMETIS
    ierr = metis_nodeND_f(nvtxs, xadj_ptr, adjncy_ptr, vwgt_ptr, options_ptr, perm_ptr, iperm_ptr)
#endif
  end function monolis_metis_nodeND

end module mod_monolis_metis
