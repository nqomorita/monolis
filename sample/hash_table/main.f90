program main
  use mod_monolis_prm
  use mod_monolis_hash
  implicit none
  integer(kind=kint) :: i, j, k, in, hash, val, ans
  character :: key*27
  logical :: is_exist, is_pushed

  write(*,"(a)")"* monolis_hash_init"
  call monolis_hash_init()

  do i = iachar("a"), iachar("k")
    do j = iachar("a"), iachar("k")
      do k = iachar("a"), iachar("k")
        val = i*1000000 + j*1000 + k
        key = trim(achar(i))//trim(achar(j))//trim(achar(k))
        call monolis_hash_push(key, val, is_pushed, is_exist)
        if(.not. is_pushed) stop "** monolis error 1"
      enddo
    enddo
  enddo

  in = 1
  do i = iachar("a"), iachar("k")
    do j = iachar("a"), iachar("k")
      do k = iachar("a"), iachar("k")
        key = trim(achar(i))//trim(achar(j))//trim(achar(k))
        ans = i*1000000 + j*1000 + k
        call monolis_hash_get(key, val, is_exist)
        if(.not. is_exist) stop "** monolis error 2"
        write(*,"(a,i12,a,a,a,i12,a,i12)")"id: ", in, ", key:", trim(key), ", val: ", val, ", ans:", ans
        if(val /= ans) stop "** monolis error 3"
        in = in + 1
      enddo
    enddo
  enddo

  write(*,"(a)")"* monolis_hash_finalize"
  call monolis_hash_finalize()
end program main
