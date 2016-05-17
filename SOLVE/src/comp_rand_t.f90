subroutine comp_rand_t (Amp, epsil, dt, taper_debut, taper_fin, duree, nstep2, spectre, nstep, Ffilt)


use fft_util

implicit none

integer, intent(IN) :: nstep, nstep2
doubleprecision, intent(IN) :: epsil, dt, taper_debut, taper_fin, duree
doubleprecision, dimension(0:2), intent(IN) :: Amp
doubleprecision, dimension(1:nstep2), intent(IN) :: spectre
doubleprecision, dimension(0:nstep-1,0:2), intent(INOUT) :: Ffilt

integer :: n, comp
doubleprecision :: hasard, wt, t0,t1,t2,t3, xi,yi, tanphi
complex*16, dimension(:,:), allocatable :: F


! Allocation
allocate (F(1:nstep2,0:2))

! Tirage des signaux aleatoires pour le GLL considere
do comp = 0,2
    if (Amp(comp)>epsil) then
        do n = 1,nstep2
            call random_number(hasard)
            hasard = hasard - 0.5
            F(n,comp) = cmplx(hasard,0.d0)
        enddo
    else
        F(:,comp) = cmplx(0.d0,0.d0)
    endif
enddo

! Filtrage
do comp = 0,2
    if (Amp(comp)>epsil) then
        call dfour1(F(:,comp),nstep2,-1)
        do n = 1,nstep2
            tanphi = aimag(F(n,comp))/real(F(n,comp))
            xi = spectre(n)/sqrt(1+tanphi**2)
            if (real(F(n,comp))<0) xi = -xi
            yi = xi*tanphi
            F(n,comp) = cmplx(xi,yi)
        enddo
        call dfour1(F(:,comp),nstep2,1)
    endif
enddo

! Taper (on attenue le debut et la fin du signal genere aleatoirement)
t0 = 0; t1 = taper_debut; t2 = duree-taper_fin; t3 = duree
do n = 0,nstep-1
    call wtcoef(n*dt,t0,t1,t2,t3,wt)
    Ffilt(n,:) = Amp(:) * wt * real(F(n+1,:))
enddo

! Desallocation
deallocate (F)


end subroutine comp_rand_t
