program pip_arg

character (len=10) :: pip_arg_ref_c, num_ref_c, size_c, dpl_c
real :: pip_arg_ref, num_ref, size, dpl, num, pi

pi=acos(-1d0)

call getarg(1,pip_arg_ref_c)
call getarg(2,num_ref_c)
call getarg(3,size_c)
call getarg(4,dpl_c)

read(pip_arg_ref_c,*) pip_arg_ref
read(num_ref_c,*) num_ref
read(size_c,*) size
read(dpl_c,*) dpl

num=4./3.*pi*(size/pi/2.*dpl)**3

print *,ceiling(pip_arg_ref*(num/num_ref)**0.33333333)

end program
