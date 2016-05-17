subroutine Define_Arrays (Tdomain, rg)


use sdomains
use angles
use tensor_util
use read_model
use module_A3d

implicit none

include 'mpif.h'

type (domain), intent (INOUT), target :: Tdomain
integer, intent(IN) :: rg

logical, parameter :: savemodel = .false.   ! retourne Vs aux GLL interieurs aux elements du chunk
                                            ! reel en coord sph (only if curve==T)
logical, parameter :: savemodel2 = .false.   ! retourne Vs aux GLL interieurs aux faces exterieures
                                             ! du chunk de ref en coord cart (only if curve==T and PML==T)
integer :: n, mat, ngllx,nglly,ngllz, ngll1,ngll2, ngll, ngllPML, ngllocean, i,j,k, n_elem, nf,ne,nv, &
           idef, code, shift, I_give_to, I_take_from, n_rings, ii,jj, meshtype, nb_interfaces
integer, parameter :: etiquette=100, epsil=1, ecart=1000
integer, dimension(mpi_status_size) :: statut
doubleprecision :: vp,vs,rho,Qmu, ri,rj,rk, dx, x,y,z, r,theta,phi, A,C,L,M,F, Gc,Gs,Hc,Hs,Bc,Bs, &
                   lambda,mu, xa,ya,za, xi,eta, Jac_surf, det, theta_ref,phi_ref, dt, lat,long, &
                   rho_anom,vp_anom,vs_anom,Qmu_anom, &
                   A_anom,C_anom,F_anom,L_anom,M_anom,Gc_anom,Gs_anom,Hc_anom,Hs_anom,Bc_anom,Bs_anom
doubleprecision, external :: pow
doubleprecision, parameter :: rhowater=1500.d0, s2=1.4142135623730950488
doubleprecision, dimension (0:1,0:1) :: LocInvGrad_surf
doubleprecision, dimension (0:2,0:2) :: Rot
doubleprecision, dimension (1:6,1:6) :: Cij
doubleprecision, dimension (:), allocatable :: timestep, rad
doubleprecision, dimension (:,:), allocatable :: coord
doubleprecision, dimension (:,:,:), allocatable :: xix,xiy,xiz, etax,etay,etaz, zetax,zetay,zetaz, Jac, Rsph, &
                                                   Rlam,Rmu,RKmod, Whei, LocMassMat, wx,wy,wz, Id
character*60 :: modelfile
character*20, parameter :: myfmt = "(4f10.3)"


!!! Attribute elastic properties !!!

meshtype = Tdomain%mesh
if (meshtype==1 .or. meshtype==2 .or. meshtype==3 .or. meshtype==6) then
   allocate (rad(0:10))
   rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
   rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
   rad(8)=6291000; rad(9)=6346600; rad(10)=6356000
else if (meshtype==5) then
   allocate (rad(0:5))
   rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5701000
   rad(4)=5971000; rad(5)=6311000
endif

if (savemodel) then
   write (modelfile,"(a,I3.3)") "model_",rg
   open (12,file=trim(modelfile))
endif
if (savemodel2) then
   write (modelfile,"(a,I3.3)") "vs_0_",rg
   open (20,file=trim(modelfile))
   write (modelfile,"(a,I3.3)") "vs_1_",rg
   open (21,file=trim(modelfile))
   write (modelfile,"(a,I3.3)") "vs_2_",rg
   open (22,file=trim(modelfile))
   write (modelfile,"(a,I3.3)") "vs_3_",rg
   open (23,file=trim(modelfile))
   write (modelfile,"(a,I3.3)") "vs_4_",rg
   open (24,file=trim(modelfile))
   write (modelfile,"(a,I3.3)") "vs_5_",rg
   open (25,file=trim(modelfile))
endif

do n = 0,Tdomain%n_elem-1
   
   mat = Tdomain%specel(n)%mat_index
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   if (Tdomain%ellipticity) then
       allocate (Rsph(0:ngllx-1,0:nglly-1,0:ngllz-1))
       call r_spheric(ngllx,nglly,ngllz,Tdomain,mat,n,Rsph)
   endif
   if (Tdomain%specel(n)%ocean==.true.) then
       allocate (Tdomain%specel(n)%hocean(0:ngllx-1,0:nglly-1))
       call find_verticale(Tdomain,n)
   endif

   do k = 0,ngllz-1
    do j = 0,nglly-1
     do i = 0,ngllx-1
        ! Taking the cartesian coordinates of the GLL
        idef = Tdomain%specel(n)%Iglobnum(i,j,k)
        x = Tdomain%Globcoord(0,idef)
        y = Tdomain%Globcoord(1,idef)
        z = Tdomain%Globcoord(2,idef)
        if (Tdomain%aniso)   call cart2sph(x,y,z,r,theta_ref,phi_ref)
        if (Tdomain%curve) then
           ! Coordinates in the real chunk
           xa = x;   ya = y;   za = z
           Rot = Tdomain%rot
           x = Rot(0,0)*xa + Rot(0,1)*ya + Rot(0,2)*za
           y = Rot(1,0)*xa + Rot(1,1)*ya + Rot(1,2)*za
           z = Rot(2,0)*xa + Rot(2,1)*ya + Rot(2,2)*za
           ! Convert the cartesian coordinates into spherical coordinates
           call cart2sph(x,y,z,r,theta,phi)
           if (Tdomain%ellipticity)   r = Rsph(i,j,k)
           if (Tdomain%specel(n)%ocean .and. k==ngllz-1) then
              Tdomain%specel(n)%hocean(i,j) = Rterre - r
              if (Tdomain%specel(n)%hocean(i,j)<0)   Tdomain%specel(n)%hocean(i,j) = 0
           endif
           ! Pay attention to the spherical interfaces
           ! L'ajustage maillage / modele de vitesse au niveau des interfaces spheriques 
           ! est pris en charge ici pour les maillages PREM.
           ! Pour tout autre maillage il est necessaire d'ajouter quelques lignes de code.
           if (meshtype==1) then
              nb_interfaces = 11
           else if (meshtype==2) then
              nb_interfaces = 8
           else if (meshtype==3) then
              nb_interfaces = 7
           else if (meshtype==5) then
              nb_interfaces = 6
           else if (meshtype==6) then
              nb_interfaces = 10
           endif
           if (meshtype==1 .or. meshtype==2 .or. meshtype==3 .or. meshtype==5 .or. meshtype==6) then
              if (k==0 .or. k==ngllz-1) then
                 interface : do idef = 0,nb_interfaces-1
                    dx = abs(r-rad(idef))
                    if (dx<ecart) then
                       r = rad(idef)
                       if (k==0)   r = r + epsil
                       if (k==ngllz-1)   r = r - epsil
                       exit interface
                    endif
                 enddo interface
              endif
           endif
        endif
        if (Tdomain%aniso) then
           if (Tdomain%curve) then
              call get_value_aniso (r,theta,phi,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Qmu,Tdomain%specel(n)%moho_position)
              !!! Les 3 lignes ci-dessous sont a decommenter ssi on veut desallouer le modele CUB
!              if (n==Tdomain%n_elem-1 .and. i==ngllx-1 .and. j==nglly-1 .and. k==ngllz-1) then
!                 call deload
!              endif
              if (savemodel) then
                 if (i/=0 .and. j/=0 .and. k/=0 .and. i/=ngllx-1 .and. j/=nglly-1 .and. k/=ngllz-1) then
                    lat = 90.d0 - rad2deg(theta)
                    long = rad2deg(phi);   if (long>180) long=long-360.d0
                    vs = dsqrt((2*L+M)/(3*rho)) ! cf Panning & Romanowicz 2006
                    write (12,*) lat, long, r/1000.d0, vs
                 endif
              endif
              if (savemodel2) then
                 xa = xa/1000.d0;   ya = ya/1000.d0;   za = za/1000.d0
                 vs = dsqrt((2*L+M)/(3*rho)) ! cf Panning & Romanowicz 2006
                 if (Tdomain%specel(n)%PML) then
                    if (Tdomain%sSubdomain(mat)%Px) then
                     if (Tdomain%sSubdomain(mat)%Left) then
                        if (i==0 .and. j/=0 .and. j/=nglly-1 .and. k/=0 .and. k/=ngllz-1)   write (24,myfmt) xa,ya,za, vs
                     else
                        if (i==ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k/=0 .and. k/=ngllz-1)   write (22,myfmt) xa,ya,za, vs
                     endif
                    endif
                    if (Tdomain%sSubdomain(mat)%Py) then
                     if (Tdomain%sSubdomain(mat)%Forward) then
                        if (i/=0 .and. i/=ngllx-1 .and. j==0 .and. k/=0 .and. k/=ngllz-1)   write (21,myfmt) xa,ya,za, vs
                     else
                        if (i/=0 .and. i/=ngllx-1 .and. j==nglly-1 .and. k/=0 .and. k/=ngllz-1)   write (23,myfmt) xa,ya,za, vs
                     endif
                    endif
                    if (Tdomain%sSubdomain(mat)%Pz) then
                     if (Tdomain%sSubdomain(mat)%Down) then
                        if (i/=0 .and. i/=ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k==0)   write (20,myfmt) xa,ya,za, vs
                     else
                        if (i/=0 .and. i/=ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k==ngllz-1)   write (25,myfmt) xa,ya,za, vs
                     endif
                    endif
                 else if (Tdomain%specel(n)%moho_position==1 .or. Tdomain%specel(n)%topo_position==1) then
                    if (i/=0 .and. i/=ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k==ngllz-1)   write (25,myfmt) xa,ya,za, vs
                 endif
              endif
              if (Tdomain%adjoint .and. (.not. Tdomain%specel(n)%PML)) then
                 Tdomain%specel(n)%save_TIparam(1,i,j,k) = A
                 Tdomain%specel(n)%save_TIparam(2,i,j,k) = C
                 Tdomain%specel(n)%save_TIparam(3,i,j,k) = F
                 Tdomain%specel(n)%save_TIparam(4,i,j,k) = L
                 Tdomain%specel(n)%save_TIparam(5,i,j,k) = M
!                 call get_value_aniso_anom (r,theta,phi,rho_anom,A_anom,C_anom,F_anom,L_anom,M_anom, &
!                                            Gc_anom,Gs_anom,Hc_anom,Hs_anom,Bc_anom,Bs_anom,Qmu_anom, &
!                                            Tdomain%specel(n)%moho_position)
                 Tdomain%specel(n)%anomaly(i,j,k) = M_anom/L_anom
              endif
           else
              call get_value_aniso (x,y,z,rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Qmu,Tdomain%specel(n)%moho_position)
           endif
           ! Ecriture du Cij en notation de Mandel et en coordonnees spheriques
           Cij(:,:) = 0.
           Cij(1,1) = C
           Cij(2,2) = A+Bc
           Cij(3,3) = A-Bc
           Cij(1,2) = F+Hc
           Cij(1,3) = F-Hc
           Cij(2,3) = A-2.d0*M
           Cij(1,4) = -s2*Hs
           Cij(2,4) = -s2*Bs/2.d0
           Cij(3,4) = -s2*Bs/2.d0
           Cij(4,4) = 2.d0*M 
           Cij(5,5) = 2.d0*(L-Gc)
           Cij(6,6) = 2.d0*(L+Gc)
           Cij(5,6) = -2.d0*Gs
           do ii = 2,6
              do jj = 1,ii-1
                 Cij(ii,jj) = Cij(jj,ii)
              enddo
           enddo
           ! Si l'elem est PML on considere un milieu isotrope moyen
           if (Tdomain%specel(n)%PML) then
              Tdomain%specel(n)%Lambda(i,j,k) = lambda_from_Cij(Cij)
              Tdomain%specel(n)%Mu(i,j,k) = mu_from_Cij(Cij)
           else ! Et sinon:
              if (Tdomain%n_sls>0) then
                 lambda = lambda_from_Cij(Cij)
                 mu = mu_from_Cij(Cij)
              endif
              ! Expression de Cij en cartesien
              call c_4tensor(Cij,theta_ref,phi_ref)
              ! Sauvegarde des 21 coeffs
              idef = 0
              do ii = 1,6
                 do jj = ii,6
                    Tdomain%specel(n)%Cij(idef,i,j,k) = Cij(ii,jj)
                    idef = idef + 1
                 enddo
              enddo
              if (Tdomain%n_sls>0) then
                 Tdomain%specel(n)%Cij( 0,i,j,k) = Tdomain%specel(n)%Cij( 0,i,j,k) - lambda-2.d0*mu
                 Tdomain%specel(n)%Cij( 6,i,j,k) = Tdomain%specel(n)%Cij( 6,i,j,k) - lambda-2.d0*mu
                 Tdomain%specel(n)%Cij(11,i,j,k) = Tdomain%specel(n)%Cij(11,i,j,k) - lambda-2.d0*mu
                 Tdomain%specel(n)%Cij( 1,i,j,k) = Tdomain%specel(n)%Cij( 1,i,j,k) - lambda
                 Tdomain%specel(n)%Cij( 2,i,j,k) = Tdomain%specel(n)%Cij( 2,i,j,k) - lambda
                 Tdomain%specel(n)%Cij( 7,i,j,k) = Tdomain%specel(n)%Cij( 7,i,j,k) - lambda
                 Tdomain%specel(n)%Cij(15,i,j,k) = Tdomain%specel(n)%Cij(15,i,j,k) - 2.d0*mu
                 Tdomain%specel(n)%Cij(18,i,j,k) = Tdomain%specel(n)%Cij(18,i,j,k) - 2.d0*mu
                 Tdomain%specel(n)%Cij(20,i,j,k) = Tdomain%specel(n)%Cij(20,i,j,k) - 2.d0*mu
                 Tdomain%specel(n)%Lambda(i,j,k) = lambda
                 Tdomain%specel(n)%Mu(i,j,k) = mu
              endif
           endif
        else
           if (Tdomain%curve) then
              call get_value (r,theta,phi,rho,vp,vs,Qmu,Tdomain%specel(n)%moho_position)
              if (savemodel) then
                 if (i/=0 .and. j/=0 .and. k/=0 .and. i/=ngllx-1 .and. j/=nglly-1 .and. k/=ngllz-1) then
                    lat = 90.d0 - rad2deg(theta)
                    long = rad2deg(phi);   if (long>180) long=long-360.d0
                    write (12,*) lat, long, r/1000.d0, vs
                 endif
              endif
              if (savemodel2) then
                 xa = xa/1000.d0;   ya = ya/1000.d0;   za = za/1000.d0
                 if (Tdomain%specel(n)%PML) then
                    if (Tdomain%sSubdomain(mat)%Px) then
                     if (Tdomain%sSubdomain(mat)%Left) then
                        if (i==0 .and. j/=0 .and. j/=nglly-1 .and. k/=0 .and. k/=ngllz-1)   write (24,myfmt) xa,ya,za, vs
                     else
                        if (i==ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k/=0 .and. k/=ngllz-1)   write (22,myfmt) xa,ya,za, vs
                     endif
                    endif
                    if (Tdomain%sSubdomain(mat)%Py) then
                     if (Tdomain%sSubdomain(mat)%Forward) then
                        if (i/=0 .and. i/=ngllx-1 .and. j==0 .and. k/=0 .and. k/=ngllz-1)   write (21,myfmt) xa,ya,za, vs
                     else
                        if (i/=0 .and. i/=ngllx-1 .and. j==nglly-1 .and. k/=0 .and. k/=ngllz-1)   write (23,myfmt) xa,ya,za, vs
                     endif
                    endif
                    if (Tdomain%sSubdomain(mat)%Pz) then
                     if (Tdomain%sSubdomain(mat)%Down) then
                        if (i/=0 .and. i/=ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k==0)   write (20,myfmt) xa,ya,za, vs
                     else
                        if (i/=0 .and. i/=ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k==ngllz-1)   write (25,myfmt) xa,ya,za, vs
                     endif
                    endif
                 else if (Tdomain%specel(n)%moho_position==1 .or. Tdomain%specel(n)%topo_position==1) then
                    if (i/=0 .and. i/=ngllx-1 .and. j/=0 .and. j/=nglly-1 .and. k==ngllz-1)   write (25,myfmt) xa,ya,za, vs
                 endif
              endif
              if (Tdomain%adjoint .and. (.not. Tdomain%specel(n)%PML)) then
!                 call get_value_anom (r,theta,phi,rho_anom,vp_anom,vs_anom,Qmu_anom,Tdomain%specel(n)%moho_position)
!                 Tdomain%specel(n)%anomaly(i,j,k) = vs_anom
              endif
           else
              call get_value (x,y,z,rho,vp,vs,Qmu,Tdomain%specel(n)%moho_position)
           endif
           Tdomain%specel(n)%Lambda(i,j,k) = (vp**2 - 2*vs**2 ) * rho
           Tdomain%specel(n)%Mu(i,j,k) = vs**2 * rho
        endif
        Tdomain%specel(n)%Density(i,j,k) = rho
        if ((.not. Tdomain%specel(n)%PML) .and. (Tdomain%n_sls>0))   Tdomain%specel(n)%Q(i,j,k) = Qmu
     enddo
    enddo
   enddo

   if (Tdomain%ellipticity)   deallocate (Rsph)

enddo

if (Tdomain%ellipticity)   deallocate(Tdomain%Rsph_Nodes)

if (savemodel) then
   close (12)
   call MPI_BARRIER(MPI_COMM_WORLD, code)
   if (rg==0) then
      call system("cat model_* > model.out")
      call system("rm model_*")
   endif
endif
if (savemodel2) then
   close (20); close (21); close (22); close (23); close (24); close (25);
   call MPI_BARRIER(MPI_COMM_WORLD, code)
   if (rg==0) then
      call system("cat vs_0_* > vs_0.out")
      call system("cat vs_1_* > vs_1.out")
      call system("cat vs_2_* > vs_2.out")
      call system("cat vs_3_* > vs_3.out")
      call system("cat vs_4_* > vs_4.out")
      call system("cat vs_5_* > vs_5.out")
      call system("rm vs_?_*")
   endif
endif


!!! Calculating the time step following the Courant criteria !!!

call compute_dt (Tdomain,dt)
allocate (timestep(0:Tdomain%n_proc-1))
call mpi_allgather(dt,1,mpi_double_precision,timestep,1,mpi_double_precision,mpi_comm_world,code)
do i = 0,Tdomain%n_proc-1
   if (timestep(i) <= dt)   dt = timestep(i)
enddo
deallocate (timestep)
Tdomain%sTimeParam%ntime = int (Tdomain%sTimeParam%Duration/dt) + 1
Tdomain%sTimeParam%dt = Tdomain%sTimeParam%Duration/(Tdomain%sTimeParam%ntime-1)
print *, "THE NUMBER OF ITERATION IS", Tdomain%sTimeParam%ntime
print *, "THE TIME STEP IS", Tdomain%sTimeParam%dt


!!! And now we can compute the source signal and allocate the traces !!!

call def_timefunc (Tdomain,rg)

if (Tdomain%save_trace) then
   do n = 0,Tdomain%n_receivers-1
      if (rg == Tdomain%sReceiver(n)%proc) then
         allocate (Tdomain%sReceiver(n)%StoreTrace (0:Tdomain%sTimeParam%ntime-1, 0:2))
         Tdomain%sReceiver(n)%StoreTrace = 0
      endif
   enddo
endif


!!! Compute Mass Matrix !!!

do n = 0,Tdomain%n_elem-1

   mat = Tdomain%specel(n)%mat_index
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   allocate (Jac (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (xix (0:ngllx-1,0:nglly-1,0:ngllz-1)) 
   allocate (xiy (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (xiz (0:ngllx-1,0:nglly-1,0:ngllz-1))    
  
   allocate (etax (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (etay (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (etaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (zetax (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (zetay (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (zetaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (Whei (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (RKmod (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (Rlam (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (Rmu (0:ngllx-1,0:nglly-1,0:ngllz-1))

   do k = 0,ngllz -1 
       do j = 0,nglly-1 
           do i = 0,ngllx-1
               Whei (i,j,k) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j) &
                              * Tdomain%sSubdomain(mat)%GLLwz(k)
               if (Tdomain%specel(n)%PML==.false.) then
                  Tdomain%specel(n)%wgtx(i) = Tdomain%sSubdomain(mat)%GLLwx(i)
                  Tdomain%specel(n)%wgty(j) = Tdomain%sSubdomain(mat)%GLLwy(j)
                  Tdomain%specel(n)%wgtz(k) = Tdomain%sSubdomain(mat)%GLLwz(k)
               endif
           enddo
       enddo
   enddo
      
   xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)
   xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
   xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)
      
   etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
   etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
   etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)

   zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
   zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
   zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)

   Jac  = Tdomain%specel(n)%Jacob

   Tdomain%specel(n)%MassMat = Whei*Tdomain%specel(n)%Density*Jac

   if (Tdomain%specel(n)%ocean) then
       if (Tdomain%n_nodes/=27) then
           print *,"THE OCEAN OPTION REQUIRES 27 CTRL_PTS FOR NOW."
           stop
       endif
       allocate (Tdomain%specel(n)%Mocean(0:ngllx-1,0:nglly-1))
       allocate (coord(0:Tdomain%n_nodes-1,0:2))
       do i = 0,Tdomain%n_nodes-1
           j = Tdomain%specel(n)%Control_Nodes(i)
           coord(i,0:2) = Tdomain%Coord_Nodes(0:2,j)
       enddo
       do j = 0,nglly-1
        eta = Tdomain%sSubdomain(mat)%GLLcy(j)
        do i = 0,ngllx-1
            xi = Tdomain%sSubdomain(mat)%GLLcx(i)
            LocInvGrad_surf(:,:) = 0
            do k = 4,7
                do ii = 0,1
                 do jj = 0,1
                     LocInvGrad_surf(ii,jj) = LocInvGrad_surf(ii,jj) + &
                                              Comp_derivshapefunc(k,xi,eta,1.d0,ii) * coord(k,jj)
                 enddo
                enddo
            enddo
            do k = 16,19
                do ii = 0,1
                 do jj = 0,1
                     LocInvGrad_surf(ii,jj) = LocInvGrad_surf(ii,jj) + &
                                              Comp_derivshapefunc(k,xi,eta,1.d0,ii) * coord(k,jj)
                 enddo
                enddo
            enddo
            do ii = 0,1
             do jj = 0,1
                 LocInvGrad_surf(ii,jj) = LocInvGrad_surf(ii,jj) + &
                                          Comp_derivshapefunc(25,xi,eta,1.d0,ii) * coord(25,jj)
             enddo
            enddo
            Jac_surf = LocInvgrad_surf(0,0) * LocInvgrad_surf(1,1) &
                       - LocInvgrad_surf(1,0) * LocInvgrad_surf(0,1)
            Tdomain%specel(n)%Mocean(i,j) = rhowater * Tdomain%specel(n)%hocean(i,j) * &
                                            Tdomain%sSubdomain(mat)%GLLwx(i) * &
                                            Tdomain%sSubdomain(mat)%GLLwy(j) * &
                                            Jac_surf
        enddo
       enddo
       deallocate (coord,Tdomain%specel(n)%hocean)
   endif
      
   if (.not. Tdomain%specel(n)%PML) then

!     Rlam = Tdomain%specel(n)%Lambda
!     Rmu  = Tdomain%specel(n)%Mu
!     RKmod = Rlam + 2.* Rmu

!     Tdomain%specel(n)%Acoeff(:,:,:,0) = -Whei*(RKmod*xix**2+Rmu*(xiy**2+xiz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,1) = -Whei*(RKmod*xix*etax+Rmu*(xiy*etay+xiz*etaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,2) = -Whei*(rKmod*xix*zetax+Rmu*(xiy*zetay+xiz*zetaz))*Jac 
!     Tdomain%specel(n)%Acoeff(:,:,:,3) = -Whei*(Rlam+Rmu)*xix*xiy*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,4) = -Whei*(Rlam*xix*etay+Rmu*xiy*etax)*Jac  
!     Tdomain%specel(n)%Acoeff(:,:,:,5) = -Whei*(Rlam*xix*zetay+Rmu*xiy*zetax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,6) = -Whei*(Rlam+Rmu)*xix*xiz*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,7) = -Whei*(Rlam*xix*etaz+Rmu*xiz*etax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,8) = -Whei*(Rlam*xix*zetaz+rmu*xiz*zetax)*Jac

!     Tdomain%specel(n)%Acoeff(:,:,:,9) = -Whei*(RKmod*etax**2+Rmu* (etay**2+etaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,10) = -Whei*(RKmod*etax*zetax+Rmu* (etay*zetay+etaz*zetaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,11) = -Whei*(Rlam*etax*xiy+Rmu*etay*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,12) = -Whei*(Rlam+Rmu)*etay*etax*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,13) = -Whei*(Rlam*etax*zetay+Rmu*etay*zetax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,14) = -Whei*(Rlam*etax*xiz+Rmu*etaz*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,15) = -Whei*(Rlam+Rmu)*etaz*etax*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,16) = -Whei*(Rlam*etax*zetaz+Rmu*etaz*zetax)*Jac

!     Tdomain%specel(n)%Acoeff(:,:,:,17) = -Whei*(RKmod*zetax**2+Rmu*  (zetay**2+zetaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,18) = -Whei*(Rlam*zetax*xiy+Rmu*zetay*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,19) = -Whei*(Rlam*zetax*etay+Rmu*zetay*etax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,20) = -Whei*(Rlam+Rmu)*zetax*zetay*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,21) = -Whei*(Rlam*zetax*xiz+Rmu*zetaz*xix)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,22) = -Whei*(Rlam*zetax*etaz+Rmu*zetaz*etax)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,23) = -Whei*(Rlam+Rmu)*zetax*zetaz*Jac

!     Tdomain%specel(n)%Acoeff(:,:,:,24) = -Whei*(RKmod*xiy**2+Rmu* (xix**2+xiz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,25) = -Whei*(RKmod*xiy*etay+Rmu* (xix*etax+xiz*etaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,26) = -Whei*(RKmod*xiy*zetay+Rmu*  (xix*zetax+xiz*zetaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei*(Rlam+Rmu)*xiy*xiz*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei*(Rlam*etaz*xiy+Rmu*etay*xiz)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei*(Rlam*zetaz*xiy+Rmu*zetay*xiz)*Jac

!     Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei*(RKmod*etay**2+Rmu* (etax**2+etaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei*(RKmod*zetay*etay+Rmu* (zetax*etax+zetaz*etaz))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei*(Rlam*etay*xiz+Rmu*etaz*xiy)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei*(Rlam+Rmu)*etay*etaz*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei*(Rlam*zetaz*etay+Rmu*zetay*etaz)*Jac

!     Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei*(RKmod*zetay**2+Rmu* (zetax**2+zetaz**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,36) = -Whei*(Rlam*xiz*zetay+Rmu*xiy*zetaz)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,37) = -Whei*(Rlam*zetay*etaz+Rmu*zetaz*etay)*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,38) = -Whei*(Rlam+Rmu)*zetay*zetaz*Jac

!     Tdomain%specel(n)%Acoeff(:,:,:,39) = -Whei*(RKmod*xiz**2+Rmu*  (xix**2+xiy**2))*Jac 
!     Tdomain%specel(n)%Acoeff(:,:,:,40) = -Whei*(RKmod*xiz*etaz+Rmu* (xix*etax+xiy*etay))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,41) = -Whei*(RKmod*xiz*zetaz+Rmu* (xix*zetax+xiy*zetay))*Jac
                                       
!     Tdomain%specel(n)%Acoeff(:,:,:,42) = -Whei*(RKmod*etaz**2+Rmu* (etax**2+etay**2))*Jac
!     Tdomain%specel(n)%Acoeff(:,:,:,43) = -Whei*(RKmod*zetaz*etaz+Rmu* (zetax*etax+zetay*etay))*Jac
                                       
!     Tdomain%specel(n)%Acoeff(:,:,:,44) = -Whei*(RKmod*zetaz**2+Rmu* (zetax**2+zetay**2))*Jac

!     deallocate (Tdomain%specel(n)%InvGrad)
                                       
   else

     RLam = Tdomain%specel(n)%Lambda
     RMu = Tdomain%specel(n)%Mu
     RKmod = RLam + 2*RMu

     Tdomain%specel(n)%Acoeff(:,:,:,0) = RKmod *xix
     Tdomain%specel(n)%Acoeff(:,:,:,1) = RKmod *etax
     Tdomain%specel(n)%Acoeff(:,:,:,2) = RKmod *zetax

     Tdomain%specel(n)%Acoeff(:,:,:,3) = RLam *xiy
     Tdomain%specel(n)%Acoeff(:,:,:,4) = RLam *etay
     Tdomain%specel(n)%Acoeff(:,:,:,5) = RLam *zetay

     Tdomain%specel(n)%Acoeff(:,:,:,6) = RLam *xiz
     Tdomain%specel(n)%Acoeff(:,:,:,7) = RLam *etaz
     Tdomain%specel(n)%Acoeff(:,:,:,8) = RLam *zetaz

     Tdomain%specel(n)%Acoeff(:,:,:,9) = RLam *xix
     Tdomain%specel(n)%Acoeff(:,:,:,10) = RLam *etax
     Tdomain%specel(n)%Acoeff(:,:,:,11) = RLam *zetax

     Tdomain%specel(n)%Acoeff(:,:,:,12) = RKmod *xiy
     Tdomain%specel(n)%Acoeff(:,:,:,13) = RKmod *etay
     Tdomain%specel(n)%Acoeff(:,:,:,14) = RKmod *zetay

     Tdomain%specel(n)%Acoeff(:,:,:,15) = RKmod *xiz
     Tdomain%specel(n)%Acoeff(:,:,:,16) = RKmod *etaz
     Tdomain%specel(n)%Acoeff(:,:,:,17) = RKmod *zetaz

     Tdomain%specel(n)%Acoeff(:,:,:,18) = RMu *xix
     Tdomain%specel(n)%Acoeff(:,:,:,19) = RMu *etax
     Tdomain%specel(n)%Acoeff(:,:,:,20) = RMu *zetax

     Tdomain%specel(n)%Acoeff(:,:,:,21) = RMu *xiy
     Tdomain%specel(n)%Acoeff(:,:,:,22) = RMu *etay
     Tdomain%specel(n)%Acoeff(:,:,:,23) = RMu *zetay

     Tdomain%specel(n)%Acoeff(:,:,:,24) = RMu *xiz
     Tdomain%specel(n)%Acoeff(:,:,:,25) = RMu *etaz
     Tdomain%specel(n)%Acoeff(:,:,:,26) = RMu *zetaz

     Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei * xix * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei * xiy * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei * xiz * Jac

     Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei * etax * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei * etay * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei * etaz * Jac

     Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei * zetax * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei * zetay * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei * zetaz * Jac

     allocate (wx (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (wy (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (wz (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (Id (0:ngllx-1,0:nglly-1,0:ngllz-1))

     if (Tdomain%sSubDomain(mat)%Px) then
         idef = Tdomain%specel(n)%Iglobnum(0,0,0)
         dx = Tdomain%GlobCoord(0,idef)
         idef = Tdomain%specel(n)%Iglobnum(ngllx-1,0,0)
         dx = abs(Tdomain%GlobCoord(0,idef) - dx) 
         if (Tdomain%sSubDomain(mat)%Left) then
            do i = 0,ngllx-1
               ri = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcx(ngllx-1-i)) * float(ngllx-1)
               vp = Rkmod(i,0,0) / Tdomain%specel(n)%Density(i,0,0)
               vp = sqrt(vp)
               wx(i,0:nglly-1,0:ngllz-1) = pow(ri, vp, ngllx-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                               Tdomain%sSubdomain(mat)%npow)
            enddo
         else 
            do i = 0,ngllx-1
               ri = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcx(i)) * float(ngllx-1)
               vp = Rkmod(i,0,0) / Tdomain%specel(n)%Density(i,0,0)
               vp = sqrt(vp)
               wx(i,0:nglly-1,0:ngllz-1) = pow(ri, vp, ngllx-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                               Tdomain%sSubdomain(mat)%npow)
            enddo
         endif
     else 
         wx = 0.
     endif  

     if (Tdomain%sSubDomain(mat)%Py) then
         idef = Tdomain%specel(n)%Iglobnum(0,0,0)
         dx = Tdomain%GlobCoord(1,idef)
         idef = Tdomain%specel(n)%Iglobnum(0,nglly-1,0)
         dx = abs(Tdomain%GlobCoord (1,idef) - dx) 
         if (Tdomain%sSubDomain(mat)%Forward) then
            do j = 0,nglly-1
               rj = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcy(nglly-1-j)) * float(nglly-1)
               vp = Rkmod(0,j,0) / Tdomain%specel(n)%Density(0,j,0)
               vp = sqrt(vp)
               wy(0:ngllx-1,j,0:ngllz-1) = pow(rj, vp, nglly-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                                Tdomain%sSubdomain(mat)%npow)
            enddo
         else 
            do j = 0,nglly-1
               rj = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcy(j)) * float(nglly-1)
               vp = Rkmod(0,j,0) / Tdomain%specel(n)%Density(0,j,0)
               vp = sqrt(vp)
               wy(0:ngllx-1,j,0:ngllz-1) = pow(rj, vp, nglly-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                                Tdomain%sSubdomain(mat)%npow)
            enddo
         endif
     else 
         wy = 0.
     endif  

     if (Tdomain%sSubDomain(mat)%Pz) then
         idef = Tdomain%specel(n)%Iglobnum(0,0,0)
         dx = Tdomain%GlobCoord(2,idef)
         idef = Tdomain%specel(n)%Iglobnum(0,0,ngllz-1)
         dx = abs(Tdomain%GlobCoord(2,idef) - dx)
         if (Tdomain%sSubDomain(mat)%Down) then
            do k = 0,ngllz-1
               rk = 0.5 * (1 + Tdomain%sSubdomain(mat)%GLLcz(ngllz-1-k)) * float(ngllz-1)
               vp = Rkmod(0,0,k) / Tdomain%specel(n)%Density(0,0,k)
               vp = sqrt(vp)
               wz(0:ngllx-1,0:nglly-1,k) = pow(rk, vp, ngllz-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                               Tdomain%sSubdomain(mat)%npow)
            enddo
         else 
            do k = 0,ngllz-1
               rk = 0.5 * (1 + Tdomain%sSubdomain(mat)%GLLcz(k)) * float(ngllz-1)
               vp = Rkmod(0,0,k) / Tdomain%specel(n)%Density(0,0,k)
               vp = sqrt(vp)
               wz(0:ngllx-1,0:nglly-1,k) = pow(rk, vp, ngllz-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                              Tdomain%sSubdomain(mat)%npow)
            enddo
         endif
     else 
         wz = 0.
     endif

     Id = 1

     Tdomain%specel(n)%DumpSx(:,:,:,1) = Id + 0.5 * Tdomain%sTimeParam%dt * wx
     Tdomain%specel(n)%DumpSx(:,:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,:,1) 
     Tdomain%specel(n)%DumpSx (:,:,:,0) = (Id - Tdomain%sTimeParam%dt * 0.5 * wx) * &
                                          Tdomain%specel(n)%DumpSx(:,:,:,1)

     Tdomain%specel(n)%DumpSy(:,:,:,1) = Id + 0.5 * Tdomain%sTimeParam%dt * wy
     Tdomain%specel(n)%DumpSy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSy (:,:,:,1) 
     Tdomain%specel(n)%DumpSy (:,:,:,0) = (Id - Tdomain%sTimeParam%dt * 0.5 * wy) * &
                                          Tdomain%specel(n)%DumpSy(:,:,:,1)

     Tdomain%specel(n)%DumpSz(:,:,:,1) = Id + 0.5 * Tdomain%sTimeParam%dt * wz
     Tdomain%specel(n)%DumpSz(:,:,:,1)  = 1./ Tdomain%specel(n)%DumpSz (:,:,:,1)
     Tdomain%specel(n)%DumpSz (:,:,:,0) = (Id - Tdomain%sTimeParam%dt * 0.5 * wz) * &
                                          Tdomain%specel(n)%DumpSz(:,:,:,1)

     Tdomain%specel(n)%DumpMass(:,:,:,0) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                                           Tdomain%sTimeParam%dt * wx * Jac
     Tdomain%specel(n)%DumpMass(:,:,:,1) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                                           Tdomain%sTimeParam%dt * wy * Jac
     Tdomain%specel(n)%DumpMass(:,:,:,2) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                                           Tdomain%sTimeParam%dt * wz * Jac

     if (Tdomain%curve) then
         call find_normales(Tdomain,n)
     endif

     deallocate (wx,wy,wz,Id)
     deallocate (Tdomain%specel(n)%InvGrad)
   endif

   deallocate (Jac, xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz, Whei, RKmod, Rmu, Rlam)  
enddo


!!! Communications !!!

do n = 0,Tdomain%n_elem-1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        select case (i)
         case (0)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,1:ngll2-2,0)
         case (1)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,0,1:ngll2-2)
         case (2)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(ngllx-1,1:ngll1-2,1:ngll2-2)
         case (3)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,nglly-1,1:ngll2-2)
         case (4)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(0,1:ngll1-2,1:ngll2-2)
         case (5)
            Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                                             Tdomain%specel(n)%MassMat(1:ngll1-2,1:ngll2-2,ngllz-1)
            if (Tdomain%sFace(nf)%ocean) then
                Tdomain%sFace(nf)%Mocean(:,:) = Tdomain%specel(n)%Mocean(1:ngll1-2,1:ngll2-2)
            endif
        end select
        if (Tdomain%sFace(nf)%PML) then
            select case (i)
             case (0)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,1:ngll2-2,0,:)
             case (1)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,0,1:ngll2-2,:)
             case (2)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,1:ngll1-2,1:ngll2-2,:)
             case (3)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,nglly-1,1:ngll2-2,:)
             case (4)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,1:ngll1-2,1:ngll2-2,:)
             case (5)
                Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll1-2,1:ngll2-2,ngllz-1,:)
            end select
        endif
    enddo
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        ngll = Tdomain%sEdge(ne)%ngll
        select case (i)
         case (0)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,0,0)
         case (1)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,1:ngll-2,0)
         case (2)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,nglly-1,0)
         case (3)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,1:ngll-2,0)
         case (4)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,0,1:ngll-2)
         case (5)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,0,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(1:ngll-2,0)
         case (6)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,0,1:ngll-2)
         case (7)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,1:ngll-2)
         case (8)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(ngllx-1,1:ngll-2,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(ngllx-1,1:ngll-2)
         case (9)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(1:ngll-2,nglly-1,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(1:ngll-2,nglly-1)
         case (10)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,nglly-1,1:ngll-2)
         case (11)
            Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                                           Tdomain%specel(n)%MassMat(0,1:ngll-2,ngllz-1)
            if (Tdomain%sEdge(ne)%ocean) Tdomain%sEdge(ne)%Mocean(:) = Tdomain%sEdge(ne)%Mocean(:) + &
                                         Tdomain%specel(n)%Mocean(0,1:ngll-2)
        end select
        if (Tdomain%sEdge(ne)%PML) then
            select case (i)
             case (0)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,0,0,:)
             case (1)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,1:ngll-2,0,:)
             case (2)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,nglly-1,0,:)
             case (3)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,1:ngll-2,0,:)
             case (4)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,0,1:ngll-2,:)
             case (5)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,0,ngllz-1,:)
             case (6)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,0,1:ngll-2,:)
             case (7)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,1:ngll-2,:)
             case (8)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(ngllx-1,1:ngll-2,ngllz-1,:)
             case (9)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(1:ngll-2,nglly-1,ngllz-1,:)
             case (10)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,nglly-1,1:ngll-2,:)
             case (11)
                Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                                                  Tdomain%specel(n)%DumpMass(0,1:ngll-2,ngllz-1,:)
            end select
        endif
    enddo
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        select case (i)
         case (0)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,0,0)
         case (1)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,0,0)
         case (2)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,0)
         case (3)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,nglly-1,0)
         case (4)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,0,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                           Tdomain%specel(n)%Mocean(0,0)
         case (5)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,0,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                            Tdomain%specel(n)%Mocean(ngllx-1,0)
         case (6)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                            Tdomain%specel(n)%Mocean(ngllx-1,nglly-1)
         case (7)
            Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                          Tdomain%specel(n)%MassMat(0,nglly-1,ngllz-1)
            if (Tdomain%sVertex(nv)%ocean) Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + &
                                           Tdomain%specel(n)%Mocean(0,nglly-1)
        end select
        if (Tdomain%sVertex(nv)%PML) then
            select case (i)
             case (0)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,0,0,:)
             case (1)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,0,0,:)
             case (2)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,0,:)
             case (3)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,nglly-1,0,:)
             case (4)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,0,ngllz-1,:)
             case (5)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,0,ngllz-1,:)
             case (6)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,ngllz-1,:)
             case (7)
                Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,nglly-1,ngllz-1,:)
            end select
        endif
    enddo
    if (Tdomain%specel(n)%ocean==.true.)   deallocate (Tdomain%specel(n)%Mocean)
enddo


!!! Invert Mass Matrix expression !!!

do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   allocate (LocMassMat(1:ngllx-2,1:nglly-2,1:ngllz-2))
   LocMassMat(:,:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2)

   if (Tdomain%specel(n)%PML) then
      Tdomain%specel(n)%DumpVx (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,0)
      Tdomain%specel(n)%DumpVx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVx (:,:,:,1) 
      Tdomain%specel(n)%DumpVx (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,0)
      Tdomain%specel(n)%DumpVx (:,:,:,0) = Tdomain%specel(n)%DumpVx (:,:,:,0) * Tdomain%specel(n)%DumpVx (:,:,:,1)

      Tdomain%specel(n)%DumpVy (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,1)
      Tdomain%specel(n)%DumpVy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVy (:,:,:,1) 
      Tdomain%specel(n)%DumpVy (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,1)
      Tdomain%specel(n)%DumpVy (:,:,:,0) = Tdomain%specel(n)%DumpVy (:,:,:,0) * Tdomain%specel(n)%DumpVy (:,:,:,1)

      Tdomain%specel(n)%DumpVz (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,2)
      Tdomain%specel(n)%DumpVz (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVz (:,:,:,1) 
      Tdomain%specel(n)%DumpVz (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,2)
      Tdomain%specel(n)%DumpVz (:,:,:,0) = Tdomain%specel(n)%DumpVz (:,:,:,0) * Tdomain%specel(n)%DumpVz (:,:,:,1)

      deallocate (Tdomain%specel(n)%DumpMass)
   endif

   LocMassmat = 1./ LocMassMat
   deallocate (Tdomain%specel(n)%MassMat) 
   allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2) )
   Tdomain%specel(n)%MassMat = LocMassMat
   deallocate (LocMassMat)
enddo


!!! MPI communications !!!

do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    ngllocean = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        do j = 1,Tdomain%sFace(nf)%ngll2-2
            do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sComm(n)%Give(ngll) = Tdomain%sFace(nf)%MassMat(k,j)
                ngll = ngll + 1
            enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2)
                    ngllPML = ngllPML + 1
                enddo
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sComm(n)%Give(ngll) = Tdomain%sEdge(ne)%MassMat(j)
            ngll = ngll + 1
        enddo
        if (Tdomain%sEdge(ne)%ocean) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%Giveocean(ngllocean) = Tdomain%sEdge(ne)%Mocean(j)
                ngllocean = ngllocean + 1
            enddo
        endif
        if (Tdomain%sEdge(ne)%PML) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2)
                ngllPML = ngllPML + 1
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        Tdomain%sComm(n)%Give(ngll) = Tdomain%svertex(nv)%MassMat
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%ocean) then
            Tdomain%sComm(n)%Giveocean(ngllocean) = Tdomain%sVertex(nv)%Mocean
            ngllocean = ngllocean + 1
        endif
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sVertex(nv)%DumpMass(0:2)
            ngllPML = ngllPML + 1
        endif
    enddo
enddo

n = Tdomain%n_proc
do shift = 1,n-1
    I_give_to = rg + shift
    if (I_give_to > n-1)   I_give_to = I_give_to - n
    I_take_from = rg - shift
    if (I_take_from < 0)   I_take_from = I_take_from + n
    if (mod(n,shift)==0 .and. shift/=1) then
        n_rings = shift
    else if (mod(n,n-shift)==0 .and. shift/=n-1) then
        n_rings = n-shift
    else if (mod(n,2)==0 .and. mod(shift,2)==0) then
        n_rings = 2
    else
        n_rings = 1
    endif
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngll>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngll>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngll>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngll>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllocean>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%Giveocean, Tdomain%sComm(I_give_to)%ngllocean, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngllocean>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%Takeocean, Tdomain%sComm(I_take_from)%ngllocean, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllocean>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%Takeocean, Tdomain%sComm(I_take_from)%ngllocean, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllocean>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%Giveocean, Tdomain%sComm(I_give_to)%ngllocean, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllPML>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngllPML>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
enddo

do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    ngllocean = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        do j = 1,Tdomain%sFace(nf)%ngll2-2
            do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sFace(nf)%MassMat(k,j) = Tdomain%sFace(nf)%MassMat(k,j) + Tdomain%sComm(n)%Take(ngll)
                ngll = ngll + 1
            enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%DumpMass(k,j,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                    ngllPML = ngllPML + 1
                enddo
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + Tdomain%sComm(n)%Take(ngll)
            ngll = ngll + 1
        enddo
        if (Tdomain%sEdge(ne)%ocean) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%Mocean(j) = Tdomain%sEdge(ne)%Mocean(j) + Tdomain%sComm(n)%Takeocean(ngllocean)
                ngllocean = ngllocean + 1
            enddo
        endif
        if (Tdomain%sEdge(ne)%PML) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%DumpMass(j,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                ngllPML = ngllPML + 1
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv = Tdomain%sComm(n)%vertices(i)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + Tdomain%sComm(n)%Take(ngll)
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%ocean) then
            Tdomain%sVertex(nv)%Mocean = Tdomain%sVertex(nv)%Mocean + Tdomain%sComm(n)%Takeocean(ngllocean)
            ngllocean = ngllocean + 1
        endif
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%sVertex(nv)%DumpMass(0:2) = Tdomain%sVertex(nv)%DumpMass(0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
            ngllPML = ngllPML + 1
        endif
    enddo
    ! A make "OPT = -O0 -C" compilation doesn't like these following deallocate. I don't know why !?!
    if (Tdomain%sComm(n)%ngll>0) then
        deallocate (Tdomain%sComm(n)%Give)
        deallocate (Tdomain%sComm(n)%Take)
    endif
    if (Tdomain%sComm(n)%ngllocean>0) then
        deallocate (Tdomain%sComm(n)%Giveocean)
        deallocate (Tdomain%sComm(n)%Takeocean)
    endif
    if (Tdomain%sComm(n)%ngllPML>0) then
        deallocate (Tdomain%sComm(n)%GivePML)
        deallocate (Tdomain%sComm(n)%TakePML)
    endif
enddo


!!! Dumping factor for PML !!!

do nf = 0,Tdomain%n_face-1
    if (Tdomain%sFace(nf)%PML) then
        Tdomain%sFace(nf)%DumpVx(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,0)
        Tdomain%sFace(nf)%DumpVx(:,:,1) = 1./Tdomain%sFace(nf)%DumpVx(:,:,1)
        Tdomain%sFace(nf)%DumpVx(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,0)
        Tdomain%sFace(nf)%DumpVx(:,:,0) = Tdomain%sFace(nf)%DumpVx(:,:,0) * Tdomain%sFace(nf)%DumpVx(:,:,1)     

        Tdomain%sFace(nf)%DumpVy(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,1) = 1./Tdomain%sFace(nf)%DumpVy (:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,0) = Tdomain%sFace(nf)%DumpVy(:,:,0) * Tdomain%sFace(nf)%DumpVy(:,:,1)         

        Tdomain%sFace(nf)%DumpVz(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,2)
        Tdomain%sFace(nf)%DumpVz(:,:,1) = 1./Tdomain%sFace(nf)%DumpVz(:,:,1)
        Tdomain%sFace(nf)%DumpVz(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,2)
        Tdomain%sFace(nf)%DumpVz(:,:,0) = Tdomain%sFace(nf)%DumpVz(:,:,0) * Tdomain%sFace(nf)%DumpVz(:,:,1)      
        deallocate (Tdomain%sFace(nf)%DumpMass)
    endif
    if (Tdomain%sFace(nf)%ocean) then
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        do i = 1,ngll1-2
         do j = 1,ngll2-2
             do ii = 0,2
              do jj = 0,2
                  Tdomain%sFace(nf)%M33ocean(i,j,ii,jj) = Tdomain%sFace(nf)%Mocean(i,j) * &
                                                          Tdomain%sFace(nf)%verticale(i,j,ii,jj)
                  if (ii==jj)   Tdomain%sFace(nf)%M33ocean(i,j,ii,jj) = Tdomain%sFace(nf)%M33ocean(i,j,ii,jj) + &
                                                                        Tdomain%sFace(nf)%MassMat(i,j)
              enddo
             enddo
             call invert_3d (Tdomain%sFace(nf)%M33ocean(i,j,0:2,0:2),det)
         enddo
        enddo
        deallocate (Tdomain%sFace(nf)%Mocean, Tdomain%sFace(nf)%verticale)
    endif
    Tdomain%sFace(nf)%MassMat = 1./ Tdomain%sFace(nf)%MassMat
enddo

do ne = 0,Tdomain%n_edge-1
    if (Tdomain%sEdge(ne)%PML) then
        Tdomain%sEdge(ne)%DumpVx(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,0)
        Tdomain%sEdge(ne)%DumpVx(:,1) = 1./Tdomain%sEdge(ne)%DumpVx(:,1)
        Tdomain%sEdge(ne)%DumpVx(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,0)
        Tdomain%sEdge(ne)%DumpVx(:,0) = Tdomain%sEdge(ne)%DumpVx(:,0) * Tdomain%sEdge(ne)%DumpVx(:,1)

        Tdomain%sEdge(ne)%DumpVy(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,1) = 1./Tdomain%sEdge(ne)%DumpVy(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,0) = Tdomain%sEdge(ne)%DumpVy(:,0) * Tdomain%sEdge(ne)%DumpVy(:,1)

        Tdomain%sEdge(ne)%DumpVz(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,2)
        Tdomain%sEdge(ne)%DumpVz(:,1) = 1./Tdomain%sEdge(ne)%DumpVz(:,1)
        Tdomain%sEdge(ne)%DumpVz(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,2)
        Tdomain%sEdge(ne)%DumpVz(:,0) = Tdomain%sEdge(ne)%DumpVz(:,0) * Tdomain%sEdge(ne)%DumpVz(:,1)
        deallocate (Tdomain%sEdge(ne)%DumpMass)
    endif
    if (Tdomain%sEdge(ne)%ocean) then
        ngll = Tdomain%sEdge(ne)%ngll
        do i = 1,ngll-2
            do ii = 0,2
             do jj = 0,2
                 Tdomain%sEdge(ne)%M33ocean(i,ii,jj) = Tdomain%sEdge(ne)%Mocean(i) * &
                                                       Tdomain%sEdge(ne)%verticale(i,ii,jj)
                 if (ii==jj)   Tdomain%sEdge(ne)%M33ocean(i,ii,jj) = Tdomain%sEdge(ne)%M33ocean(i,ii,jj) + &
                                                                     Tdomain%sEdge(ne)%MassMat(i)
             enddo
            enddo
            call invert_3d (Tdomain%sEdge(ne)%M33ocean(i,0:2,0:2),det)
        enddo
        deallocate (Tdomain%sEdge(ne)%Mocean, Tdomain%sEdge(ne)%verticale)
    endif
    Tdomain%sEdge(ne)%MassMat = 1./ Tdomain%sEdge(ne)%MassMat
enddo

do nv = 0,Tdomain%n_vertex-1
    if (Tdomain%sVertex(nv)%PML) then
        Tdomain%sVertex(nv)%DumpVx(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(0)
        Tdomain%sVertex(nv)%DumpVx(1) = 1./Tdomain%sVertex(nv)%DumpVx(1)
        Tdomain%sVertex(nv)%DumpVx(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(0)
        Tdomain%sVertex(nv)%DumpVx(0) = Tdomain%sVertex(nv)%DumpVx(0) * Tdomain%sVertex(nv)%DumpVx(1)

        Tdomain%sVertex(nv)%DumpVy(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(1)
        Tdomain%sVertex(nv)%DumpVy(1) = 1./Tdomain%sVertex(nv)%DumpVy(1)
        Tdomain%sVertex(nv)%DumpVy(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(1)
        Tdomain%sVertex(nv)%DumpVy(0) = Tdomain%sVertex(nv)%DumpVy(0) * Tdomain%sVertex(nv)%DumpVy(1)

        Tdomain%sVertex(nv)%DumpVz(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(2)
        Tdomain%sVertex(nv)%DumpVz(1) = 1./Tdomain%sVertex(nv)%DumpVz(1)
        Tdomain%sVertex(nv)%DumpVz(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(2)
        Tdomain%sVertex(nv)%DumpVz(0) = Tdomain%sVertex(nv)%DumpVz(0) * Tdomain%sVertex(nv)%DumpVz(1)
        deallocate (Tdomain%sVertex(nv)%DumpMass)
    endif
    if (Tdomain%sVertex(nv)%ocean) then
        do ii = 0,2
         do jj = 0,2
             Tdomain%sVertex(nv)%M33ocean(ii,jj) = Tdomain%sVertex(nv)%Mocean * &
                                                   Tdomain%sVertex(nv)%verticale(ii,jj)
             if (ii==jj)   Tdomain%sVertex(nv)%M33ocean(ii,jj) = Tdomain%sVertex(nv)%M33ocean(ii,jj) + &
                                                                 Tdomain%sVertex(nv)%MassMat
         enddo
        enddo
        call invert_3d (Tdomain%sVertex(nv)%M33ocean(0:2,0:2),det)
        deallocate (Tdomain%sVertex(nv)%verticale)
    endif
    Tdomain%sVertex(nv)%MassMat = 1./ Tdomain%sVertex(nv)%MassMat
enddo


return
end subroutine

