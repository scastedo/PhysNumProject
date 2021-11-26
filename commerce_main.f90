module param
  implicit none
  integer, parameter :: nat = 127, T0 = 1e3  ! nombre de villes/points/...
  real(8), dimension(nat,2) :: xvec, xvec_new  ! positions des villes/points (2d)
  end module param


program main
  use param
  implicit none
  real(8) :: distance, diff, temperature_function  ! fonctions
  integer, allocatable :: seed(:)

  integer :: m, istep
  integer :: nstep = 10000000  ! nombre d'itérations
  real(8), parameter :: kB = 0.08617  ! constante de Boltzmann [meV/K]
  real(8) :: T  ! température fictive [K]
  integer :: compteur_ta  ! compteur de tirages acceptés

  real(8) :: D, dD
  real(8) :: s1, s2, s  ! nombre aléatoire
  integer :: k1, k2  ! indice des points à échanger
    

  call random_seed(size=m); allocate(seed(m)) ; seed = 14119265 ; call random_seed(put=seed) !fix random seed
  	
! initialiser xvec
  call initialize_xvec()

  ! écrire dans un fichier les positions initiales (xvec)
  open(1, file='pos_init_beer.res')  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(2, file="dist_beer.res") !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(3, file='pos_fin_beer.res')  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do m = 1, nat
     write(1,*) xvec(m,1), xvec(m,2)
  enddo
  close(1)



  ! calculer et écrire à l'écran la distance initiale
  D = distance()
  write(*,*) 'Distance initiale =', D, 'm'

  
  write(2,*) '0', D
  
  ! initialiser le compteur de tirages acceptés
  compteur_ta = 0
  
  ! boucle Monte-Carlo
  do istep = 1, nstep

     T = temperature_function(istep, nstep)
     
     ! tirage des points au hasard
     call random_number(s1)
     call random_number(s2)
     k1 = floor(nat*s1)+1
     k2= floor(nat*s2)+1
     ! tirage des points à échanger
     xvec_new = xvec
     xvec_new(k1,1) = xvec(k2,1)
     xvec_new(k1,2) = xvec(k2,2)
     xvec_new(k2,1) = xvec(k1,1)
     xvec_new(k2,2) = xvec(k1,2)
     
     ! calcul de la différence de distance dD induite par l'échange
     dD = diff(k1, k2)
     
     ! application du critère de Metropolis
     if (dD.lt.0) then
        xvec = xvec_new
        D = D + dD 
        compteur_ta = compteur_ta + 1
     
     else
        call random_number(s)
        if (s.le.exp(-dD/(kB*T))) then
           xvec = xvec_new
           D = D + dD 
           compteur_ta = compteur_ta + 1
        endif        
        
     endif
     
     ! écrire (tous les 1000 pas) la nouvelle distance dans un fichier
     if (mod(istep,1000).eq.0) then
        write(2,*) istep, D
     endif
    
  enddo
  
  close(2)
  
  ! écrire à l'écran la distance finale et le taux d'acceptation des tirages
  write(*,*) "distance finale =", D, "m"
  write(*,*) "taux d'acceptation des tirages =", real(compteur_ta)/nstep
  
  ! écrire dans un fichier les positions finales
  do m = 1, nat
     write(3,*) xvec(m,1), xvec(m,2)
  enddo
  close(3)
  
end program main


subroutine initialize_xvec()
  use param
  implicit none
  real(8) :: b1,b2
  integer :: b, c1, c2,i
  real(8), dimension(1,2) :: temp
  
  open(4, file = 'bier127.txt')
  do i = 1, nat
     read(4,*) xvec(i,1), xvec(i,2)
  enddo
  
  do b = 0 , (nat*5)
     call random_number(b1)
     call random_number(b2)
     c1 = floor(nat*b1+1)
     c2= floor(nat*b2+1)
     temp(1,1) = xvec(c1,1)
     temp(1,2) = xvec(c1,2)
     xvec(c1,1) = xvec(c2,1)
     xvec(c1,2) = xvec(c2,2)
     xvec(c2,1) = temp(1,1)
     xvec(c2,2) = temp(1,2)
  enddo
  close(4)
end subroutine initialize_xvec


real(8) function diff(k1, k2)
   use param
   implicit none
   integer :: k1, k2, i
   real(8) :: distance, distance_new
   diff = 0.0
   distance = 0.0
   distance_new = 0.0

   do i = 1, nat-1
	      distance = distance + sqrt((xvec(i,1)-xvec(i+1,1))**2 + (xvec(i,2)-xvec(i+1,2))**2)
        enddo
        distance = distance +  sqrt((xvec(nat,1)-xvec(1,1))**2 + (xvec(nat,2)-xvec(1,2))**2)


   do i = 1, nat-1
	      distance_new = distance_new + sqrt((xvec_new(i,1)-xvec_new(i+1,1))**2 + (xvec_new(i,2)-xvec_new(i+1,2))**2)
        enddo
        distance_new = distance_new +  sqrt((xvec_new(nat,1)-xvec_new(1,1))**2 + (xvec_new(nat,2)-xvec_new(1,2))**2)

   diff = distance_new - distance


!!$   diff = diff + sqrt((xvec_new(k1,1) - xvec_new((k1-1) +floor(nat*1.0/(nat*1.0+(k1-1)))*nat,1))**2 + (xvec_new(k1,2) -&
!!$        xvec_new((k1-1) +floor(nat*1.0/(nat*1.0+(k1-1)))*nat,2))**2)
!!$   diff = diff + sqrt((xvec_new(mod(k1+1,nat),1) - xvec_new(k1,1))**2 + (xvec_new(mod(k1+1,nat),2) - xvec_new(k1,2))**2)
!!$   diff = diff + sqrt((xvec_new(k2,1) - xvec_new((k2-1) +floor(nat*1.0/(nat*1.0+(k2-1)))*nat,1))**2 + (xvec_new(k2,2) -&
!!$        xvec_new((k2-1) +floor(nat*1.0/(nat*1.0+(k2-1)))*nat,2))**2)
!!$   diff = diff + sqrt((xvec_new(mod(k2+1,nat),1) - xvec_new(k2,1))**2 + (xvec_new(mod(k2+1,nat),2) - xvec_new(k2,2))**2)
!!$
!!$   diff = diff - sqrt((xvec(k1,1) - xvec((k1-1) +floor(nat*1.0/(nat*1.0+(k1-1)))*nat,1))**2 + (xvec(k1,2) - xvec((k1-1)&
!!$        +floor(nat*1.0/(nat*1.0+(k1-1)))*nat,2))**2)
!!$   diff = diff - sqrt((xvec(mod(k1+1,nat),1) - xvec(k1,1))**2 + (xvec(mod(k1+1,nat),2) - xvec(k1,2))**2)
!!$   diff = diff - sqrt((xvec(k2,1) - xvec((k2-1) +floor(nat*1.0/(nat*1.0+(k2-1)))*nat,1))**2 + (xvec(k2,2) - xvec((k2-1)&
!!$        +floor(nat*1.0/(nat*1.0+(k2-1)))*nat,2))**2)
!!$   diff = diff - sqrt((xvec(mod(k2+1,nat),1) - xvec(k2,1))**2 + (xvec(mod(k2+1,nat),2) - xvec(k2,2))**2)
end function diff


! (histogram)




real(8) function distance()
	use param
	implicit none
	integer :: i, j
	distance = 0.0
	do i = 1, nat-1
	      distance = distance + sqrt((xvec(i,1)-xvec(i+1,1))**2 + (xvec(i,2)-xvec(i+1,2))**2)
        enddo
        distance = distance +  sqrt((xvec(nat,1)-xvec(1,1))**2 + (xvec(nat,2)-xvec(1,2))**2)
end function distance



real(8) function temperature_function(istep, nstep)
  use param
  implicit none
  integer :: istep, nstep
  real :: linear_temp, exp_temp, const_temp
  linear_temp = T0 - istep*T0/nstep
  exp_temp   = T0*((0.8)**istep) !exponential multiplicative cooling Kirkpatrick, Gelatt and Vecchi (1983)
  const_temp = 0.5
  temperature_function = linear_temp
end function temperature_function

        

subroutine init_random_seed()
	use iso_fortran_env, only: int64
	implicit none
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid
	integer(int64) :: t

	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(newunit=un, file="/dev/urandom", access="stream", &
	form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
	   read(un) seed
	   close(un)
	else
	! Fallback to XOR:ing the current time and pid. The PID is
	! useful in case one launches multiple instances of the same
	! program in parallel.
	   call system_clock(t)
	   if (t == 0) then
	      call date_and_time(values=dt)
	      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
		  + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
		  + dt(3) * 24_int64 * 60 * 60 * 1000 &
		  + dt(5) * 60 * 60 * 1000 &
		  + dt(6) * 60 * 1000 + dt(7) * 1000 + dt(8)
	   end if
	   pid = getpid()
	   t = ieor(t, int(pid, kind(t)))
	   do i = 1, n
	      seed(i) = lcg(t)
	   end do
	end if
	call random_seed(put=seed)
	     contains
		function lcg(s)
		   integer :: lcg
		   integer(int64) :: s
		   if (s == 0) then
		      s = 104729
		   else
		      s = mod(s, 4294967296_int64)
		   end if
		   s = mod(s * 279470273_int64, 4294967291_int64)
		   lcg = int(mod(s, int(huge(0), int64)), kind(0))
		end function lcg
end subroutine init_random_seed
