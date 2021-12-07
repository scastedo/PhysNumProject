module param
   implicit none
   integer, parameter :: n_cities = 48, T0 = 3000  ! nombre de villes, température initiale
   real(8), dimension(n_cities,2) :: xvec, xvec_new  ! positions (2d) des villes dans l'ordre du chemin
end module param



program main
   use param
   implicit none
   
   real(8) :: distance, diff, temperature_function  ! fonctions

   integer :: m, istep
   integer :: nstep = 100000  ! nombre d'itérations
   real(8), parameter :: kB = 1  ! constante de Boltzmann
   real(8) :: T  ! température fictive [K]
   integer :: compteur_ta  ! compteur de tirages acceptés

   real(8) :: D, dD  ! distance, différence en distance entre deux parcours
   real(8) :: s1, s2, s  ! nombres aléatoires
   integer :: k1, k2  ! indices des points à échanger
  	
   ! initialiser xvec
   call initialize_xvec()

   ! ouvrir des fichiers pour y écrire quelques résultats plus tard
   open(1, file='pos_init_us.res')
   open(2, file='dist_us.res') 
   open(3, file='pos_fin_us.res')  

   ! écrire dans un fichier les positions initiales
   do m = 1, n_cities
      write(1,*) xvec(m,1), xvec(m,2)
   enddo
   close(1)

   ! calculer et écrire à l'écran la distance initiale
   D = distance()
   write(*,*) "distance initiale =", D, "km"
  
   write(2,*) "0", D
  
   ! initialiser le compteur de tirages acceptés
   compteur_ta = 0
  
  
   !!! BOUCLE MONTE-CARLO !!!
   do istep = 1, nstep

      T = temperature_function(istep, nstep)
     
      ! tirage des points au hasard
      call random_number(s1)
      call random_number(s2)
      k1 = floor(n_cities*s1)+1
      k2 = floor(n_cities*s2)+1
     
      ! tirage des points à échanger dans l'ordre
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
   write(*,*) "distance finale =", D, "km"
   write(*,*) "taux d'acceptation des tirages =", real(compteur_ta)/nstep
  
   ! écrire dans un fichier les positions finales
   do m = 1, n_cities
      write(3,*) xvec(m,1), xvec(m,2)
   enddo
   close(3)
  
end program main



! ------------------------------------------------------------------------------------
! --------- subroutine pour lire la carte et initialiser un chemin randomisé ---------
! ------------------------------------------------------------------------------------
subroutine initialize_xvec()
   use param
   implicit none
   real(8) :: b1, b2  ! nombres aléatoires
   integer :: b, i, c1, c2  ! integers pour boucles et nombres aléatoires
   real(8), dimension(1,2) :: temp
  
   ! lire les coordonnées des villes dans un ordre initial
   open(4, file='datasets/us48.txt')
   do i = 1, n_cities
      read(4,*) xvec(i,1), xvec(i,2)
   enddo
  
   ! randomiser l'ordre
   do b = 0 , (n_cities*5)
      call random_number(b1)
      call random_number(b2)
      c1 = floor(n_cities*b1 + 1)
      c2 = floor(n_cities*b2 + 1)
      
      temp(1,1) = xvec(c1,1)
      temp(1,2) = xvec(c1,2)
      xvec(c1,1) = xvec(c2,1)
      xvec(c1,2) = xvec(c2,2)
      xvec(c2,1) = temp(1,1)
      xvec(c2,2) = temp(1,2)
   enddo
  
   close(4)
  
end subroutine initialize_xvec



! ------------------------------------------------------------------------------------
! ---- fonction pour calculer la différence de longueur totale entre deux chemins ----
! ----        dont les villes à positions k1 et k2 ont été échangées              ----
! ------------------------------------------------------------------------------------             
real(8) function diff(k1, k2)
   use param
   implicit none
   integer :: k1, k2, i
   real(8) :: distance, distance_new
   
   diff = 0.0
   distance = 0.0
   distance_new = 0.0

   do i = 1, n_cities-1
	      distance = distance + sqrt((xvec(i,1)-xvec(i+1,1))**2 + (xvec(i,2)-xvec(i+1,2))**2)
   enddo
   distance = distance +  sqrt((xvec(n_cities,1)-xvec(1,1))**2 + (xvec(n_cities,2)-xvec(1,2))**2)


   do i = 1, n_cities-1
	      distance_new = distance_new + sqrt((xvec_new(i,1)-xvec_new(i+1,1))**2 + (xvec_new(i,2)-xvec_new(i+1,2))**2)
   enddo
   distance_new = distance_new +  sqrt((xvec_new(n_cities,1)-xvec_new(1,1))**2 + (xvec_new(n_cities,2)-xvec_new(1,2))**2)

   diff = distance_new - distance

end function diff


! ------------------------------------------------------------------------------------
! -------------- fonction pour calculer la longueur totale d'un chemin ---------------
! ------------------------------------------------------------------------------------
real(8) function distance()
	 use param
 	 implicit none
	 integer :: i
 	 
   distance = 0.0
 	 
   do i = 1, n_cities-1
	    distance = distance + sqrt((xvec(i,1)-xvec(i+1,1))**2 + (xvec(i,2)-xvec(i+1,2))**2)
   enddo
   distance = distance + sqrt((xvec(n_cities,1)-xvec(1,1))**2 + (xvec(n_cities,2)-xvec(1,2))**2)

end function distance



! ------------------------------------------------------------------------------------
! ---------- fonciton pour calculer la température au pas d'itération istep ----------
! -  selon le schéma choisi (constant, pas à pas, linéaire, exponentiel ou sigmoide) -
! ------------------------------------------------------------------------------------
real(8) function temperature_function(istep, nstep)
   use param
   implicit none
   integer :: istep, nstep
   real :: linear_temp, exp_temp, const_temp, stepwise_temp, sigmoid_temp
 
   linear_temp = T0 - istep * (T0-1)/nstep
   exp_temp   = T0 * ((0.8)**(istep/5000))  ! exponential multiplicative cooling Kirkpatrick, Gelatt and Vecchi (1983)
   const_temp = T0
   sigmoid_temp = 2500 / (0.5 + exp((istep*1.0-40000)/10000))

   ! stepwise cooling
   if (istep.lt.nstep/9) then
      stepwise_temp = 3000
   else if (istep.lt.(2*nstep/9).and.istep.ge.(nstep/9)) then
      stepwise_temp = 3000 - (istep-nstep/9)*2000/(nstep/9)
   else if (istep.lt.nstep/3.and.istep.ge.(2*nstep/9)) then
      stepwise_temp = 1000
   else if (istep.lt.(4*nstep/9).and.istep.ge.(nstep/3)) then
         stepwise_temp = 1000 - (istep-nstep/3)*950/(nstep/9)
   else if (istep.lt.(5*nstep/9).and.istep.ge.(4*nstep/9)) then 
       stepwise_temp = 50
    else if (istep.lt.(2*nstep/3).and.istep.ge.(5*nstep/9)) then
       stepwise_temp = 50 + (istep-5*nstep/9)*350/(nstep/9)
    else if (istep.lt.(7*nstep/9).and.istep.ge.(2*nstep/3)) then
       stepwise_temp = 400
    else if (istep.lt.(8*nstep/9).and.istep.ge.(7*nstep/9)) then
        stepwise_temp = 400 - (istep-7*nstep/9)*398/(nstep/9)
    else if (istep.ge.(8*nstep/9)) then
       stepwise_temp = 2
    endif
   
  
   temperature_function = sigmoid_temp  ! choisir la température/le cooling

end function temperature_function