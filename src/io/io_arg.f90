module mod_monolis_io_arg
  use mod_monolis_prm
  use mod_monolis_com
  use mod_monolis_util
  use mod_monolis_util_debug
  use mod_monolis_mesh
  use mod_monolis_stdlib
  implicit none

contains

  subroutine monolis_get_part_arg(n_domain, is_format_id, is_overlap)
    implicit none
    integer(kint) :: i, count, n, n_domain
    character :: argc1*128, argc2*128
    logical :: is_overlap, is_format_id

    call monolis_debug_header("monolis_get_part_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h")then
        write(*,"(a)")"-n {num of subdomain}: the number of subdomain"
        write(*,"(a)")"-t {N/O}: type of domain decomposition (N:non-overlapping, O:overlapping)"
        write(*,"(a)")"--with-id {Y/N}: node or elem id appears at the beginning of each line"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    n_domain = 1
    is_overlap = .true.
    is_format_id = .false.

    if(mod(count,2) /= 0) stop "* monolis partitioner input arg error"
    do i = 1, count/2
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)
      if(trim(argc1) == "-n")then
        read(argc2,*) n
        n_domain = n

      elseif(trim(argc1) == "-t")then
        if(trim(argc2) == "O") is_overlap = .true.
        if(trim(argc2) == "N") is_overlap = .false.

      elseif(trim(argc1) == "--with-id")then
        if(trim(argc2) == "Y") is_format_id = .true.
        if(trim(argc2) == "N") is_format_id = .false.

      else
        write(*,"(a)")"* monolis input arg error"
        stop
      endif
    enddo

    call monolis_debug_int("n_domain", n_domain)
  end subroutine monolis_get_part_arg

  subroutine monolis_get_nodal_graph_part_arg(fname, n_domain)
    implicit none
    integer(kint) :: i, count, n, n_domain
    character :: argc1*128, argc2*128, fname*100

    call monolis_debug_header("monolis_get_nodal_graph_part_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h")then
        write(*,"(a)")"-n {num of subdomain}"
        write(*,"(a)")"-i {input file name}"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    n_domain = 1
    fname = "graph.dat"

    if(mod(count,2) /= 0) stop "* monolis partitioner input arg error"
    do i = 1, count/2
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)
      if(trim(argc1) == "-n")then
        read(argc2,*) n
        n_domain = n

      elseif(trim(argc1) == "-i")then
        fname = trim(argc2)

      else
        write(*,"(a)")"* monolis input arg error"
        stop
      endif
    enddo

    call monolis_debug_int("n_domain", n_domain)
  end subroutine monolis_get_nodal_graph_part_arg

  subroutine monolis_get_connectivity_part_arg(fname, n_domain)
    implicit none
    integer(kint) :: i, count, n, n_domain
    character :: argc1*128, argc2*128, fname*100

    call monolis_debug_header("monolis_get_connectivity_part_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h")then
        write(*,"(a)")"-n {num of subdomain}"
        write(*,"(a)")"-i {input file name}"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    n_domain = 1
    fname = "connectivity.dat"

    if(mod(count,2) /= 0) stop "* monolis partitioner input arg error"
    do i = 1, count/2
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)
      if(trim(argc1) == "-n")then
        read(argc2,*) n
        n_domain = n

      elseif(trim(argc1) == "-i")then
        fname = trim(argc2)

      else
        write(*,"(a)")"* monolis input arg error"
        stop
      endif
    enddo

    call monolis_debug_int("n_domain", n_domain)
  end subroutine monolis_get_connectivity_part_arg

  subroutine monolis_get_part_bc_arg(n_domain, fname)
    implicit none
    integer(kint) :: i, count, n, n_domain
    character :: argc1*128, argc2*128, fname*100

    call monolis_debug_header("monolis_get_part_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h")then
        write(*,"(a)")"-n {num of subdomain}: the number of subdomain"
        write(*,"(a)")"-i {input filename}: input filename"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    n_domain = 1
    fname = "D_bc.dat"

    if(mod(count,2) /= 0) stop "* monolis partitioner input arg error"
    do i = 1, count/2
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)

      if(trim(argc1) == "-n")then
        read(argc2,*) n
        n_domain = n
      elseif(trim(argc1) == "-i")then
        fname = trim(argc2)
      else
        write(*,"(a)")"* monolis input arg error"
        stop
      endif
    enddo

    call monolis_debug_int("n_domain", n_domain)
  end subroutine monolis_get_part_bc_arg

  subroutine monolis_get_dbc_all_arg(n_block, val, fnname, fename, foname)
    implicit none
    integer(kint) :: i, j, count, n, n_block
    real(kdouble), allocatable :: val(:)
    character :: argc1*128, argc2*128, fnname*100, fename*100, foname*100

    call monolis_debug_header("monolis_get_dbc_all_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 0 .or. count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h" .or. count == 0)then
        write(*,"(a)")"usage:"
        write(*,"(a)") &
        & "./monolis_dbc_all {options} {block size} {value 1} {value 2} ... {value n}"
        write(*,"(a)")""
        write(*,"(a)")"options:"
        write(*,"(a)")"-in {input node filename}: (defualt) node.dat"
        write(*,"(a)")"-ie {input elem filename}: (defualt) elem.dat"
        write(*,"(a)")"-o  {output filename}: (defualt) D_bc.dat"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    fnname = "node.dat"
    fename = "elem.dat"
    foname = "D_bc.dat"

    j = 0
    do i = 1, count/2
      j = i
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)

      if(trim(argc1) == "-in")then
        fnname = trim(argc2)
      elseif(trim(argc1) == "-ie")then
        fename = trim(argc2)
      elseif(trim(argc1) == "-o")then
        foname = trim(argc2)
      else
        exit
      endif
    enddo

    call getarg(2*j-1, argc1)
    read(argc1,*) n_block
    allocate(val(n_block), source = 0.0d0)

    do i = 1, n_block
      call getarg(2*j-1 + i, argc1)
      read(argc1,*) val(i)
    enddo
  end subroutine monolis_get_dbc_all_arg

  subroutine monolis_get_extract_all_arg(n_block, val, fnname, fename, foname)
    implicit none
    integer(kint) :: i, j, count, n, n_block
    real(kdouble), allocatable :: val(:)
    character :: argc1*128, argc2*128, fnname*100, fename*100, foname*100

    call monolis_debug_header("monolis_get_extract_all_arg")

    call monolis_set_debug(.true.)

    count = iargc()
    if(count == 1)then
      call getarg(1, argc1)
      if(trim(argc1) == "-h")then
        write(*,"(a)")"usage:"
        write(*,"(a)") &
        & "./monolis_dbc_all {options} {block size} {value 1} {value 2} ... {value n}"
        write(*,"(a)")""
        write(*,"(a)")"options:"
        write(*,"(a)")"-in {input node filename}: (defualt) node.dat"
        write(*,"(a)")"-ie {input elem filename}: (defualt) elem.dat"
        write(*,"(a)")"-o  {output filename}: (defualt) D_bc.dat"
        write(*,"(a)")"-h: help"
        stop
      endif
    endif

    fnname = "node.dat"
    fename = "elem.dat"
    foname = "surf.dat"

    j = 0
    do i = 1, count/2
      j = i
      call getarg(2*i-1, argc1)
      call getarg(2*i  , argc2)

      if(trim(argc1) == "-in")then
        fnname = trim(argc2)
      elseif(trim(argc1) == "-ie")then
        fename = trim(argc2)
      elseif(trim(argc1) == "-o")then
        foname = trim(argc2)
      else
        exit
      endif
    enddo
  end subroutine monolis_get_extract_all_arg

end module mod_monolis_io_arg
