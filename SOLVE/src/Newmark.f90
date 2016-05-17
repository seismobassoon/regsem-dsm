subroutine Newmark (Tdomain,rg,ntime)


use sdomains
use angles

implicit none

include 'mpif.h'

type (domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: rg, ntime

logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
integer :: n, mat, elem, ngll1,ngll2,ngll3, ngllx,nglly,ngllz, nf,ne,nv, ngll, code, ntime_rand_F, &
           i,j,k, x,y,z, shift, I_give_to, I_take_from, n_rings, ngllPML, nfbis, i_simu
integer, parameter :: etiquette = 100
integer, dimension(mpi_status_size) :: statut
doubleprecision :: dt, alpha,bega,gam1, ft, cosg,sing,longdiff, u,v,w, xa,ya,za, r,theta,phi, &
                   lat_deg,phi_deg, dt_rand_F, weight, integrale
doubleprecision, dimension(0:2) :: tmp, tmp2, ft_adj
doubleprecision, dimension(:), allocatable :: get_integrale
doubleprecision, dimension (0:2,0:2) :: Rot
character*60 :: fnamef, snapfile, adjfile


! Predictor-multicorrector Newmark velocity scheme within a 
! time staggered stress-velocity formulation inside PML


dt = Tdomain%sTimeParam%dt


! Predictor Phase   

alpha = Tdomain%sTimeParam%alpha
bega = Tdomain%sTimeParam%beta / Tdomain%sTimeParam%gamma
gam1 = 1. / Tdomain%sTimeParam%gamma

do i_simu = 0,Tdomain%nb_simu-1

 do n = 0,Tdomain%n_elem-1 
     mat = Tdomain%specel(n)%mat_index
     if (.not. Tdomain%specel(n)%PML) then
         call Prediction_Elem_Veloc (Tdomain%specel(n), alpha, bega, gam1, dt, i_simu)
     else
         call get_PMLprediction_v2el (Tdomain, n, bega, dt, i_simu)
         call get_PMLprediction_e2el (Tdomain, n, bega, dt, i_simu)
         call get_PMLprediction_f2el (Tdomain, n, bega, dt, i_simu)
         if (Tdomain%curve) then
             call Prediction_Elem_PML_Veloc_curve (Tdomain%specel(n), bega, dt, &
                                                   Tdomain%sSubDomain(mat)%hTprimex, &
                                                   Tdomain%sSubDomain(mat)%hprimey, &
                                                   Tdomain%sSubDomain(mat)%hprimez, &
                                                   Tdomain%specel(n)%sSimu(i_simu))
         else
             call Prediction_Elem_PML_Veloc (Tdomain%specel(n), bega, dt, &
                                             Tdomain%sSubDomain(mat)%hTprimex, &
                                             Tdomain%sSubDomain(mat)%hprimey, &
                                             Tdomain%sSubDomain(mat)%hprimez, &
                                             Tdomain%specel(n)%sSimu(i_simu))
         endif
     endif
 enddo

 allocate (L_Face(0:Tdomain%n_face-1))
 L_Face = .true.
 allocate (L_Edge(0:Tdomain%n_edge-1))
 L_Edge = .true.
 allocate (L_Vertex(0:Tdomain%n_vertex-1))
 L_Vertex = .true.
 do n = 0,Tdomain%n_elem-1
     do i = 0,5
         nf = Tdomain%specel(n)%Near_Faces(i)
         if (L_Face(nf)) then
             L_Face(nf) = .false.
             if (.not.Tdomain%sface(nf)%PML) &
              call Prediction_Face_Veloc (Tdomain%sface(nf), alpha, bega, gam1, dt, i_simu)
         endif
     enddo
     do i = 0,11
         ne = Tdomain%specel(n)%Near_Edges(i)
         if (L_Edge(ne)) then
             L_Edge(ne) = .false.
             if (.not.Tdomain%sedge(ne)%PML) &
              call Prediction_Edge_Veloc (Tdomain%sedge(ne), alpha, bega, gam1, dt, i_simu)
         endif
     enddo
     do i = 0,7
         nv = Tdomain%specel(n)%Near_Vertices(i)
         if (L_Vertex(nv)) then
             L_Vertex(nv) = .false.
             if (.not.Tdomain%svertex(nv)%PML) &
              call Prediction_Vertex_Veloc (Tdomain%svertex(nv), alpha, bega, gam1, dt, i_simu)
         endif
     enddo
 enddo
 deallocate (L_Face,L_Edge,L_Vertex)

enddo


if (Tdomain%t_reversal_mirror==1 .and. ntime/=0)   call rec_mirror (Tdomain, ntime, rg)
if (Tdomain%t_reversal_mirror==2)   call impose_mirror (Tdomain, ntime, rg)


! Solution phase

do i_simu = 0,Tdomain%nb_simu-1

 do n = 0,Tdomain%n_elem-1
     mat = Tdomain%specel(n)%mat_index
     if (.not. Tdomain%specel(n)%PML) then
        call get_Displ_Face2Elem (Tdomain, n, i_simu)
        call get_Displ_Edge2Elem (Tdomain, n, i_simu)
        call get_Displ_Vertex2Elem (Tdomain, n, i_simu)
        call forces_int (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
                         Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                         Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                         Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls, Tdomain%aniso, Tdomain%adjoint)
!!! Si on veut utiliser ce qui suit a la place de la routine forces_int, alors il faut
!!! decommenter l'allocation de Acoeff dans "allocate_domain.f90" et la definition de
!!! Acoeff dans "define_arrays.f90".
!!! ATTENTION: la routine ci-dessous ne prend en compte ni l'attenuation ni l'anisotropie.
!        call compute_InternalForces_Elem (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
!                                          Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
!                                          Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
!                                          Tdomain%sSubDomain(mat)%hTprimez)
     else
        call compute_InternalForces_PML_Elem (Tdomain%specel(n), i_simu, Tdomain%sSubDomain(mat)%hprimex, &
                                              Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hTprimez)
     endif
 enddo

enddo

! External Forces

do n = 0, Tdomain%n_source-1
   if (rg == Tdomain%sSource(n)%proc) then
    if (Tdomain%sSource(n)%i_type_source==1 .or. Tdomain%sSource(n)%i_type_source==2) then
       elem = Tdomain%Ssource(n)%elem
       if (Tdomain%sSource(n)%i_time_function == 3) then
          ft = Tdomain%sSource(n)%timefunc(ntime)
          if (Tdomain%t_reversal_mirror==2)   ft = Tdomain%sSource(n)%timefunc(Tdomain%sTimeParam%ntime-1-ntime)
       else
          ft = CompSource(Tdomain%sSource(n),Tdomain%sTimeParam%rtime,Tdomain%t_reversal_mirror,Tdomain%sTimeParam%Duration)
          if (Tdomain%t_reversal_mirror==2)   ft = CompSource(Tdomain%sSource(n), &
                                                              Tdomain%sTimeParam%Duration-Tdomain%sTimeParam%rtime, &
                                                              Tdomain%t_reversal_mirror,Tdomain%sTimeParam%Duration)
       endif
       do x = 0,Tdomain%specel(elem)%ngllx-1
        do y = 0,Tdomain%specel(elem)%nglly-1
         do z = 0,Tdomain%specel(elem)%ngllz-1
            Tdomain%specel(elem)%sSimu(0)%Forces(x,y,z,:) = Tdomain%specel(elem)%sSimu(0)%Forces(x,y,z,:) + &
                                                            ft * Tdomain%sSource(n)%coeff(x,y,z,:)
         enddo
        enddo
       enddo
    endif
   endif
   if (Tdomain%sSource(n)%i_type_source==4) then
      dt_rand_F = Tdomain%sTimeParam%dt_rand_F
      ntime_rand_F = int(ntime*dt/dt_rand_F)
      weight = (ntime*dt - ntime_rand_F*dt_rand_F) / dt_rand_F
      do elem = 0,Tdomain%n_elem-1
         if (Tdomain%specel(elem)%random_t) then
            do x = 0,Tdomain%specel(elem)%ngllx-1
             do y = 0,Tdomain%specel(elem)%nglly-1
                z = Tdomain%specel(elem)%ngllz-1
                tmp(:) = weight * Tdomain%specel(elem)%random_coeff(x,y,ntime_rand_F+1,:) + &
                         (1.d0-weight) * Tdomain%specel(elem)%random_coeff(x,y,ntime_rand_F,:)
                Tdomain%specel(elem)%sSimu(0)%Forces(x,y,z,:) = Tdomain%specel(elem)%sSimu(0)%Forces(x,y,z,:) &
                                                                + tmp(:)
             enddo
            enddo
         endif
      enddo
   endif
enddo

! Adjoint sources

if (Tdomain%adjoint) then
   do n = 0, Tdomain%n_source_adj-1
      if (rg == Tdomain%sAdjoint(n)%proc) then
         ft_adj(0) = Tdomain%sAdjoint(n)%timefunc(0,ntime)
         ft_adj(1) = Tdomain%sAdjoint(n)%timefunc(1,ntime)
         ft_adj(2) = Tdomain%sAdjoint(n)%timefunc(2,ntime)
         elem = Tdomain%sAdjoint(n)%elem
         do x = 0,Tdomain%specel(elem)%ngllx-1
          do y = 0,Tdomain%specel(elem)%nglly-1
           do z = 0,Tdomain%specel(elem)%ngllz-1
              Tdomain%specel(elem)%sSimu(1)%Forces(x,y,z,:) = Tdomain%specel(elem)%sSimu(1)%Forces(x,y,z,:) + &
                                                              ft_adj(:) * Tdomain%sAdjoint(n)%coeff(x,y,z)
           enddo
          enddo
         enddo
      endif
   enddo
endif


! Communication of Forces

do i_simu = 0,Tdomain%nb_simu-1

 do nf = 0,Tdomain%n_face-1
    Tdomain%sFace(nf)%sSimu(i_simu)%Forces = 0
    if (Tdomain%sFace(nf)%PML) then
       Tdomain%sFace(nf)%sSimu(i_simu)%Forces1 = 0
       Tdomain%sFace(nf)%sSimu(i_simu)%Forces2 = 0
       Tdomain%sFace(nf)%sSimu(i_simu)%Forces3 = 0
    endif
 enddo
 do ne = 0,Tdomain%n_edge-1
    Tdomain%sEdge(ne)%sSimu(i_simu)%Forces = 0
    if (Tdomain%sEdge(ne)%PML) then
       Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1 = 0
       Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2 = 0
       Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3 = 0
    endif
 enddo
 do nv = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(nv)%sSimu(i_simu)%Forces = 0
    if (Tdomain%sVertex(nv)%PML) then
       Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1 = 0
       Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2 = 0
       Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3 = 0
    endif
 enddo

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
           Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0,0:2)
           if (Tdomain%sFace(nf)%PML) then
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0,0:2)
           endif
        case (1)
           Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll1-2,0,1:ngll2-2,0:2)
           if (Tdomain%sFace(nf)%PML) then
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll1-2,0,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll1-2,0,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll1-2,0,1:ngll2-2,0:2)
           endif
        case (2)
           Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                         Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
           if (Tdomain%sFace(nf)%PML) then
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,1:ngll1-2,1:ngll2-2,0:2)
           endif
        case (3)
          Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                        Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
           if (Tdomain%sFace(nf)%PML) then
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll1-2,nglly-1,1:ngll2-2,0:2)
           endif
        case (4)
           Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces(0,1:ngll1-2,1:ngll2-2,0:2)
           if (Tdomain%sFace(nf)%PML) then
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                               Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,1:ngll1-2,1:ngll2-2,0:2)
           endif
        case (5)
           Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                         Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
           if (Tdomain%sFace(nf)%PML) then
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
            Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll1-2,1:ngll2-2,ngllz-1,0:2)
           endif
       end select
    enddo
    do i = 0,11
       ne = Tdomain%specel(n)%Near_Edges(i)
       ngll = Tdomain%sEdge(ne)%ngll
       select case (i)
        case (0)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll-2,0,0,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll-2,0,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll-2,0,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll-2,0,0,0:2)
           endif
        case (1)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,1:ngll-2,0,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,1:ngll-2,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,1:ngll-2,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,1:ngll-2,0,0:2)
           endif
        case (2)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll-2,nglly-1,0,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll-2,nglly-1,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll-2,nglly-1,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll-2,nglly-1,0,0:2)
           endif
        case (3)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(0,1:ngll-2,0,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,1:ngll-2,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,1:ngll-2,0,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,1:ngll-2,0,0:2)
           endif
        case (4)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,0,1:ngll-2,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,0,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,0,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,0,1:ngll-2,0:2)
           endif
        case (5)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll-2,0,ngllz-1,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll-2,0,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll-2,0,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll-2,0,ngllz-1,0:2)
           endif
        case (6)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(0,0,1:ngll-2,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,0,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,0,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,0,1:ngll-2,0:2)
           endif
        case (7)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,nglly-1,1:ngll-2,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,nglly-1,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,nglly-1,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,nglly-1,1:ngll-2,0:2)
           endif
        case (8)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,1:ngll-2,ngllz-1,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,1:ngll-2,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,1:ngll-2,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,1:ngll-2,ngllz-1,0:2)
           endif
        case (9)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngll-2,nglly-1,ngllz-1,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(1:ngll-2,nglly-1,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(1:ngll-2,nglly-1,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(1:ngll-2,nglly-1,ngllz-1,0:2)
           endif
        case (10)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(0,nglly-1,1:ngll-2,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,nglly-1,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,nglly-1,1:ngll-2,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,nglly-1,1:ngll-2,0:2)
           endif
        case (11)
           Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(1:ngll-2,0:2) + &
                                                                  Tdomain%specel(n)%sSimu(i_simu)%Forces(0,1:ngll-2,ngllz-1,0:2)
           if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,1:ngll-2,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,1:ngll-2,ngllz-1,0:2)
            Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(1:ngll-2,0:2) + &
                                                                    Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,1:ngll-2,ngllz-1,0:2)
           endif
       end select
    enddo
    do i = 0,7
       nv = Tdomain%specel(n)%Near_Vertices(i)
       select case (i)
        case (0)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(0,0,0,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,0,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,0,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,0,0,0:2)
           endif
        case (1)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,0,0,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,0,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,0,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,0,0,0:2)
           endif
        case (2)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,nglly-1,0,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,nglly-1,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,nglly-1,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,nglly-1,0,0:2)
           endif
        case (3)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(0,nglly-1,0,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,nglly-1,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,nglly-1,0,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,nglly-1,0,0:2)
           endif
        case (4)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(0,0,ngllz-1,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,0,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,0,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,0,ngllz-1,0:2)
           endif
        case (5)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,0,ngllz-1,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,0,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,0,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,0,ngllz-1,0:2)
           endif
        case (6)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,nglly-1,ngllz-1,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(ngllx-1,nglly-1,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(ngllx-1,nglly-1,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(ngllx-1,nglly-1,ngllz-1,0:2)
           endif
        case (7)
           Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                           Tdomain%specel(n)%sSimu(i_simu)%Forces(0,nglly-1,ngllz-1,0:2)
           if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces1(0,nglly-1,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces2(0,nglly-1,ngllz-1,0:2)
            Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                             Tdomain%specel(n)%sSimu(i_simu)%Forces3(0,nglly-1,ngllz-1,0:2)
           endif
       end select
    enddo
 enddo

 ! MPI communications

 do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
       nf = Tdomain%sComm(n)%faces(i)
       do j = 1,Tdomain%sFace(nf)%ngll2-2
          do k = 1,Tdomain%sFace(nf)%ngll1-2
             Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(k,j,0:2)
             ngll = ngll + 1
          enddo
       enddo
       if (Tdomain%sFace(nf)%PML) then
          do j = 1,Tdomain%sFace(nf)%ngll2-2
             do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(k,j,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(k,j,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(k,j,0:2)
                ngllPML = ngllPML + 1
             enddo
          enddo
       endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
       ne = Tdomain%sComm(n)%edges(i)
       do j = 1,Tdomain%sEdge(ne)%ngll-2
          Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(j,0:2)
          ngll = ngll + 1
       enddo
       if (Tdomain%sEdge(ne)%PML) then
          do j = 1,Tdomain%sEdge(ne)%ngll-2
             Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(j,0:2)
             Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(j,0:2)
             Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(j,0:2)
             ngllPML = ngllPML + 1
          enddo
       endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
       nv = Tdomain%sComm(n)%vertices(i)
       Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2)
       ngll = ngll + 1
       if (Tdomain%sVertex(nv)%PML) then
          Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2)
          Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2)
          Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2)
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
           call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForces, 3*Tdomain%sComm(I_give_to)%ngll, &
                          MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
          endif
          if (Tdomain%sComm(I_take_from)%ngll>0) then
           call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForces, 3*Tdomain%sComm(I_take_from)%ngll, &
                          MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
          endif
       else
          do j = 0,n/n_rings-1
             if (rg == i + j*n_rings) then
                if (Tdomain%sComm(I_take_from)%ngll>0) then
                   call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForces, 3*Tdomain%sComm(I_take_from)%ngll, &
                                  MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                endif
                if (Tdomain%sComm(I_give_to)%ngll>0) then
                   call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForces, 3*Tdomain%sComm(I_give_to)%ngll, &
                                  MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                endif
             endif
          enddo
       endif
    enddo
    do i = 0,n_rings-1
       if (rg==i) then
          if (Tdomain%sComm(I_give_to)%ngllPML>0) then
           call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesPML, 9*Tdomain%sComm(I_give_to)%ngllPML, &
                          MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
          endif
          if (Tdomain%sComm(I_take_from)%ngllPML>0) then
           call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesPML, 9*Tdomain%sComm(I_take_from)%ngllPML, &
                          MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
          endif
       else
          do j = 0,n/n_rings-1
             if (rg == i + j*n_rings) then
                if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                 call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesPML, 9*Tdomain%sComm(I_take_from)%ngllPML, &
                                MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                endif
                if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                 call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesPML, 9*Tdomain%sComm(I_give_to)%ngllPML, &
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
    do i = 0,Tdomain%sComm(n)%nb_faces-1
       nf = Tdomain%sComm(n)%faces(i)
       do j = 1,Tdomain%sFace(nf)%ngll2-2
          do k = 1,Tdomain%sFace(nf)%ngll1-2
             Tdomain%sFace(nf)%sSimu(i_simu)%Forces(k,j,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces(k,j,0:2) + &
                                                               Tdomain%sComm(n)%TakeForces(ngll,0:2)
             ngll = ngll + 1
          enddo
       enddo
       if (Tdomain%sFace(nf)%PML) then
          do j = 1,Tdomain%sFace(nf)%ngll2-2
             do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(k,j,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces1(k,j,0:2) + &
                                                                   Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(k,j,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces2(k,j,0:2) + &
                                                                   Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(k,j,0:2) = Tdomain%sFace(nf)%sSimu(i_simu)%Forces3(k,j,0:2) + &
                                                                   Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                ngllPML = ngllPML + 1
             enddo
          enddo
       endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
       ne = Tdomain%sComm(n)%edges(i)
       do j = 1,Tdomain%sEdge(ne)%ngll-2
          Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(j,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces(j,0:2) + &
                                                          Tdomain%sComm(n)%TakeForces(ngll,0:2)
          ngll = ngll + 1
       enddo
       if (Tdomain%sEdge(ne)%PML) then
          do j = 1,Tdomain%sEdge(ne)%ngll-2
             Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(j,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces1(j,0:2) + &
                                                              Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
             Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(j,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces2(j,0:2) + &
                                                              Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
             Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(j,0:2) = Tdomain%sEdge(ne)%sSimu(i_simu)%Forces3(j,0:2) + &
                                                              Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
             ngllPML = ngllPML + 1
          enddo
       endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
       nv =  Tdomain%sComm(n)%vertices(i)
       Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces(0:2) + &
                                                       Tdomain%sComm(n)%TakeForces(ngll,0:2)
       ngll = ngll + 1
       if (Tdomain%sVertex(nv)%PML) then
          Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces1(0:2) + &
                                                           Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
          Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces2(0:2) + &
                                                           Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
          Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) = Tdomain%sVertex(nv)%sSimu(i_simu)%Forces3(0:2) + &
                                                           Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
          ngllPML = ngllPML + 1
       endif
    enddo
 enddo

enddo


! Correction phase

do i_simu = 0,Tdomain%nb_simu-1

 allocate (L_Face(0:Tdomain%n_face-1))
 L_Face = .true.
 allocate (L_Edge(0:Tdomain%n_edge-1))
 L_Edge = .true.
 allocate (L_Vertex(0:Tdomain%n_vertex-1))
 L_Vertex = .true.
 do n = 0,Tdomain%n_elem-1
    if (.not. Tdomain%specel(n)%PML) then
       call Correction_Elem_Veloc (Tdomain%specel(n), i_simu, bega, gam1, dt)
    else 
       call Correction_Elem_PML_Veloc (Tdomain%specel(n), i_simu, dt)
    endif
    do i = 0,5
       nf = Tdomain%specel(n)%Near_Faces(i)
       if (L_Face(nf)) then
          L_Face(nf) = .false.
          if (.not. Tdomain%sface(nf)%PML) then
             call Correction_Face_Veloc (Tdomain%sface(nf), i_simu, bega, gam1, dt)
          else
             call Correction_Face_PML_Veloc (Tdomain%sface(nf), i_simu, dt)
          endif
       endif
    enddo
    do i = 0,11
       ne = Tdomain%specel(n)%Near_Edges(i)
       if (L_Edge(ne)) then
          L_Edge(ne) = .false.
          if (.not. Tdomain%sedge(ne)%PML) then
             call Correction_Edge_Veloc (Tdomain%sedge(ne), i_simu, bega, gam1, dt)
          else
             call Correction_Edge_PML_Veloc (Tdomain%sedge(ne), i_simu, dt)
          endif
       endif
    enddo
    do i = 0,7
       nv = Tdomain%specel(n)%Near_Vertices(i)
       if (L_Vertex(nv)) then
          L_Vertex(nv) = .false.
          if (.not. Tdomain%svertex(nv)%PML) then
             call Correction_Vertex_Veloc (Tdomain%svertex(nv), i_simu, bega, gam1, dt)
          else
             call Correction_Vertex_PML_Veloc (Tdomain%svertex(nv), i_simu, dt)
          endif
       endif
    enddo
 enddo
 deallocate (L_Face,L_Edge,L_Vertex)

enddo


! Save Trace

if (Tdomain%save_trace) then
   do n = 0, Tdomain%n_receivers-1
      if (rg == Tdomain%sReceiver(n)%proc) then
         i = Tdomain%sReceiver(n)%elem

         ngll1 = Tdomain%specel(i)%ngllx
         ngll2 = Tdomain%specel(i)%nglly
         ngll3 = Tdomain%specel(i)%ngllz
         do x = 0,ngll1-1
          do y = 0,ngll2-1
           do z = 0,ngll3-1
              if (x==0) then
                 if (y==0) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(0)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(4)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(6)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(z,:)
                    endif
                 else if (y==ngll2-1) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(3)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(7)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(10)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(z,:)
                    endif
                 else if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(3)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(y,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(11)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(y,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(4)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%sSimu(0)%Displ(y,z,:)
                 endif
              else if (x==ngll1-1) then
                 if (y==0) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(1)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(5)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(4)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(z,:)
                    endif
                 else if (y==ngll2-1) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(2)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(6)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%sSimu(0)%Displ(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(7)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(z,:)
                    endif
                 else if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(1)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(y,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(8)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(y,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(2)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%sSimu(0)%Displ(y,z,:)
                 endif
              else if (y==0) then
                 if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(0)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(x,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(5)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(x,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(1)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%sSimu(0)%Displ(x,z,:)
                 endif
              else if (y==ngll2-1) then
                 if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(2)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(x,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(9)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%sSimu(0)%Displ(x,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(3)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%sSimu(0)%Displ(x,z,:)
                 endif
              else if (z==0) then
                 nf = Tdomain%specel(i)%Near_Faces(0)
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%sSimu(0)%Displ(x,y,:)
              else if (z==ngll3-1) then
                 nf = Tdomain%specel(i)%Near_Faces(5)
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%sSimu(0)%Displ(x,y,:)
              else
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%specel(i)%sSimu(0)%Displ(x,y,z,:)
              endif
           enddo
          enddo
         enddo
         do x = 0,ngll1-1
          do y = 0,ngll2-1
           do z = 0,ngll3-1
              Tdomain%sReceiver(n)%StoreTrace(ntime,:) = Tdomain%sReceiver(n)%StoreTrace(ntime,:) + &
                                Tdomain%sReceiver(n)%coeff(x,y,z,:) * Tdomain%sReceiver(n)%pol(x,y,z)
           enddo
          enddo
         enddo

         if (Tdomain%curve) then
          if (Tdomain%comp_rot) then
             do i = 0,2
                tmp(i) = 0
                do j = 0,2
                   tmp(i) = tmp(i) + Tdomain%sReceiver(n)%Passage(i,j) * Tdomain%sReceiver(n)%StoreTrace(ntime,j)
                enddo
             enddo
             do i = 0,2
                Tdomain%sReceiver(n)%StoreTrace(ntime,i) = tmp(i)
             enddo
          else
             do i = 0,2
                tmp(i) = 0
                do j = 0,2
                   tmp(i) = tmp(i) + Tdomain%rot(i,j) * Tdomain%sReceiver(n)%StoreTrace(ntime,j)
                enddo
             enddo
             do i = 0,2
                tmp2(i) = 0
                do j = 0,2
                   tmp2(i) = tmp2(i) + Tdomain%sReceiver(n)%Pass(i,j) * tmp(j)
                enddo
             enddo
             do i = 0,2
                Tdomain%sReceiver(n)%StoreTrace(ntime,i) = tmp2(i)
             enddo
          endif
         endif

         if (ntime==0) then
            if (Tdomain%curve) then
               if (Tdomain%comp_rot) then
                  write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_Z"
                  open (61+n*3,file=fnamef,status="replace",form="formatted")
                  write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_R"
                  open (62+n*3,file=fnamef,status="replace",form="formatted")
                  write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_T"
                  open (63+n*3,file=fnamef,status="replace",form="formatted")
               else
                  write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_Z"
                  open (61+n*3,file=fnamef,status="replace",form="formatted")
                  write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_N"
                  open (62+n*3,file=fnamef,status="replace",form="formatted")
                  write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_E"
                  open (63+n*3,file=fnamef,status="replace",form="formatted")
               endif
            else
               write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_X"
               open (61+n*3,file=fnamef,status="replace",form="formatted")
               write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_Y"
               open (62+n*3,file=fnamef,status="replace",form="formatted")
               write (fnamef,"(a,a)") trim(Tdomain%sReceiver(n)%sta_name),"_Z"
               open (63+n*3,file=fnamef,status="replace",form="formatted")
            endif
         endif

         if (Tdomain%curve .and. Tdomain%comp_rot) then
            cosg = Tdomain%sReceiver(n)%cosgamma
            sing = Tdomain%sReceiver(n)%singamma
            longdiff = Tdomain%sSource(0)%reflong - Tdomain%sReceiver(n)%reflong
            if (longdiff < 0)   longdiff = longdiff + 2*pi
            tmp(:) = Tdomain%sReceiver(n)%StoreTrace(ntime,:)
            if (longdiff < pi) then ! l'onde va d'Est en Ouest
               Tdomain%sReceiver(n)%StoreTrace(ntime,1) = -cosg*tmp(1) -sing*tmp(2)
               Tdomain%sReceiver(n)%StoreTrace(ntime,2) = sing*tmp(1) -cosg*tmp(2)
            else ! l'onde va d'Ouest en Est
               Tdomain%sReceiver(n)%StoreTrace(ntime,1) = -cosg*tmp(1) +sing*tmp(2)
               Tdomain%sReceiver(n)%StoreTrace(ntime,2) = -sing*tmp(1) -cosg*tmp(2)
            endif
         endif

         if (Tdomain%sTimeParam%samp_period<dt .or. ntime==0) then
            write (61+n*3,*) ntime*dt, Tdomain%sReceiver(n)%StoreTrace(ntime,0)
            write (62+n*3,*) ntime*dt, Tdomain%sReceiver(n)%StoreTrace(ntime,1)
            write (63+n*3,*) ntime*dt, Tdomain%sReceiver(n)%StoreTrace(ntime,2)
         else    
            if (Tdomain%sTimeParam%rtime .GT. Tdomain%sTimeParam%t_resamp) then
               weight = (Tdomain%sTimeParam%rtime - Tdomain%sTimeParam%t_resamp) / dt
               tmp(:) = (1.d0-weight)*Tdomain%sReceiver(n)%StoreTrace(ntime,:) + &
                         weight*Tdomain%sReceiver(n)%StoreTrace(ntime-1,:)
               write (61+n*3,*) Tdomain%sTimeParam%t_resamp, tmp(0)
               write (62+n*3,*) Tdomain%sTimeParam%t_resamp, tmp(1)
               write (63+n*3,*) Tdomain%sTimeParam%t_resamp, tmp(2)
            endif
         endif

         if (ntime==Tdomain%sTimeParam%ntime-1) then
            close (61+n*3)
            close (62+n*3)
            close (63+n*3)
         endif

      endif
   enddo
   if (Tdomain%sTimeParam%rtime .GT. Tdomain%sTimeParam%t_resamp)   Tdomain%sTimeParam%t_resamp = Tdomain%sTimeParam%t_resamp + &
                                                                                   Tdomain%sTimeParam%samp_period
endif


! Kernel computation (interaction between the adjoint and regular wavefields)

if (Tdomain%adjoint) then
   call kernel_rho (Tdomain,dt,ntime)
   if (Tdomain%aniso) then
      call kernel_aniso (Tdomain,dt,ntime)
   else
      call kernel_iso (Tdomain,dt,ntime)
   endif
   if (ntime==Tdomain%sTimeParam%ntime-1) then
      call save_kernel (Tdomain,rg,integrale)   ! computes and saves the vs kernel (and also the xi kernel if there is anisotropy)
      ! Petit check de la derivee partielle
!      call MPI_BARRIER(MPI_COMM_WORLD, code)
!      allocate (get_integrale(0:Tdomain%n_proc-1))
!      call mpi_gather(integrale,1,mpi_double_precision,get_integrale,1,mpi_double_precision,0,mpi_comm_world,code)
!      if (rg==0) then
!         do n = 1,Tdomain%n_proc-1
!            get_integrale(0) = get_integrale(0) + get_integrale(n)
!         enddo
!         print *, "DELTA_CHI = ", get_integrale(0)
!      endif
!      deallocate(get_integrale)
   endif
endif


! Save Snapshots

if (Tdomain%save_snapshots .and. mod(ntime,Tdomain%sTimeParam%dnsnap)==0) then
   write (snapfile,"(a,I3.3,a,I3.3)") "snapshot",Tdomain%sTimeParam%Nsnap,"_",rg
   open (12,file=trim(snapfile))
   if (Tdomain%adjoint) then
      write (adjfile,"(a,I3.3,a,I3.3)") "adjsnap_",Tdomain%sTimeParam%Nsnap,"_",rg
      open (22,file=trim(adjfile))
   endif
   do n = 0,Tdomain%n_elem-1
      if (.not. Tdomain%specel(n)%PML) then
      ngll1 = Tdomain%specel(n)%ngllx
      ngll2 = Tdomain%specel(n)%nglly
      ngll3 = Tdomain%specel(n)%ngllz
      do z = 0,ngll3-1
       do y = 0,ngll2-1
        do x = 0,ngll1-1
           i = Tdomain%specel(n)%Iglobnum(x,y,z)
           u = Tdomain%Globcoord(0,i)
           v = Tdomain%Globcoord(1,i)
           w = Tdomain%Globcoord(2,i)
           if (Tdomain%curve) then
              ! Coordinates in the real chunk
              xa = u;   ya = v;   za = w
              Rot = Tdomain%rot
              u = Rot(0,0)*xa + Rot(0,1)*ya + Rot(0,2)*za
              v = Rot(1,0)*xa + Rot(1,1)*ya + Rot(1,2)*za
              w = Rot(2,0)*xa + Rot(2,1)*ya + Rot(2,2)*za
              ! Convert cartesian coordinates into spherical coordinates
              call cart2sph(u,v,w,r,theta,phi)
              lat_deg = 90 - 180*theta/pi
              phi_deg = 180*phi/pi
              if (phi_deg>180)   phi_deg = phi_deg-360
              tmp(0) = lat_deg;   tmp(1) = phi_deg;   tmp(2) = r/1000
           else
              tmp(0) = u/1000;   tmp(1) = v/1000;   tmp(2) = w/1000
           endif
           if (x/=0 .and. y/=0 .and. z/=0 .and. x/=ngll1-1 .and. y/=ngll2-1 .and. z/=ngll3-1) then
              write (12,*) tmp(0:2)
              write (12,*) dsqrt(Tdomain%specel(n)%sSimu(0)%Displ(x,y,z,0)**2 + &
                                 Tdomain%specel(n)%sSimu(0)%Displ(x,y,z,1)**2 + &
                                 Tdomain%specel(n)%sSimu(0)%Displ(x,y,z,2)**2)
              if (Tdomain%adjoint) then
                 write (22,*) tmp(0:2)
                 write (22,*) dsqrt(Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,0)**2 + &
                                    Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,1)**2 + &
                                    Tdomain%specel(n)%sSimu(1)%Displ(x,y,z,2)**2)
              endif
           endif
        enddo
       enddo
      enddo
      endif
   enddo
   close (12)
   if (Tdomain%adjoint) then
      close (22)
   endif
endif


return
end subroutine Newmark
