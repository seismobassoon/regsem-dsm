subroutine forces_int(Elem, i_simu, hprimex, htprimex, hprimey, htprimey, hprimez, htprimez, n_solid, aniso, adj)


use sdomains

implicit none

type (Element), intent (INOUT) :: Elem
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex, htprimex
doubleprecision, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
doubleprecision, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
integer, intent(IN) :: n_solid, i_simu
logical, intent(IN) :: aniso, adj

integer :: n_z, m1,m2,m3, i,j,k, mat
doubleprecision :: epsilon_trace_over_3
doubleprecision, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
                                                                               dUy_dxi, dUy_deta, dUy_dzeta, &
                                                                               dUz_dxi, dUz_deta, dUz_dzeta, &
                                                                               DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
                                                                               Fox,Foy,Foz
doubleprecision, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
                                                  epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc


m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Elem%sSimu(i_simu)%Forces(0,0,0,0), m1, 0.d0, dUx_dxi, m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,n_z,0), m1, hprimey, m2, 0.d0, dUx_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,0,0), m1*m2, hprimez, m3, 0.d0, dUx_dzeta, m1*m2)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Elem%sSimu(i_simu)%Forces(0,0,0,1), m1, 0.d0, dUy_dxi, m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,n_z,1), m1, hprimey, m2, 0.d0, dUy_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,0,1), m1*m2, hprimez, m3, 0.d0, dUy_dzeta, m1*m2)

call DGEMM ('N', 'N', m1, m2*m3, m1, 1.d0, htprimex, m1, Elem%sSimu(i_simu)%Forces(0,0,0,2), m1, 0.d0, dUz_dxi, m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM ('N', 'N', m1, m2, m2, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,n_z,2), m1, hprimey, m2, 0.d0, dUz_deta(0,0,n_z), m1)
enddo
call DGEMM ('N', 'N', m1*m2, m3, m3, 1.d0, Elem%sSimu(i_simu)%Forces(0,0,0,2), m1*m2, hprimez, m3, 0.d0, dUz_dzeta, m1*m2)

do i = 0,m1-1
 do j = 0,m2-1
  do k = 0,m3-1

     dxx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
     dyy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)
     dzz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

     dyx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)
     dzx(i,j,k) = dUx_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUx_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUx_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

     dxy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
     dzy(i,j,k) = dUy_dxi(i,j,k)*Elem%InvGrad(i,j,k,2,0) + dUy_deta(i,j,k)*Elem%InvGrad(i,j,k,2,1) + dUy_dzeta(i,j,k)*Elem%InvGrad(i,j,k,2,2)

     dxz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,0,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,0,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,0,2)
     dyz(i,j,k) = dUz_dxi(i,j,k)*Elem%InvGrad(i,j,k,1,0) + dUz_deta(i,j,k)*Elem%InvGrad(i,j,k,1,1) + dUz_dzeta(i,j,k)*Elem%InvGrad(i,j,k,1,2)

  enddo
 enddo
enddo

if (adj) then
   do i = 0,m1-1
    do j = 0,m2-1
     do k = 0,m3-1
        Elem%sSimu(i_simu)%save_strain(i,j,k,0) = DXX(i,j,k)
        Elem%sSimu(i_simu)%save_strain(i,j,k,1) = DYY(i,j,k)
        Elem%sSimu(i_simu)%save_strain(i,j,k,2) = DZZ(i,j,k)
        Elem%sSimu(i_simu)%save_strain(i,j,k,3) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
        Elem%sSimu(i_simu)%save_strain(i,j,k,4) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
        Elem%sSimu(i_simu)%save_strain(i,j,k,5) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
     enddo
    enddo
   enddo
endif

if (n_solid>0) then
   allocate (epsilondev_xx_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_yy_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_xy_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_xz_loc(0:m1-1,0:m2-1,0:m3-1))
   allocate (epsilondev_yz_loc(0:m1-1,0:m2-1,0:m3-1))
   do i = 0,m1-1
    do j = 0,m2-1
     do k = 0,m3-1
        epsilon_trace_over_3 = 0.3333333333333333d0 * (DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k))
        epsilondev_xx_loc(i,j,k) = DXX(i,j,k) - epsilon_trace_over_3
        epsilondev_yy_loc(i,j,k) = DYY(i,j,k) - epsilon_trace_over_3
        epsilondev_xy_loc(i,j,k) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
        epsilondev_xz_loc(i,j,k) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
        epsilondev_yz_loc(i,j,k) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
     enddo
    enddo
   enddo
endif

if (aniso) then
   if (n_solid>0) then
      call calcul_forces_aniso_att(Fox,Foy,Foz, &
                                   Elem%Invgrad(:,:,:,0,0), &
                                   Elem%Invgrad(:,:,:,1,0), &
                                   Elem%Invgrad(:,:,:,2,0), &
                                   Elem%Invgrad(:,:,:,0,1), &
                                   Elem%Invgrad(:,:,:,1,1), &
                                   Elem%Invgrad(:,:,:,2,1), &
                                   Elem%Invgrad(:,:,:,0,2), &
                                   Elem%Invgrad(:,:,:,1,2), &
                                   Elem%Invgrad(:,:,:,2,2), &
                                   htprimex, htprimey, htprimez, &
                                   Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                                   DXX,DXY,DXZ, &
                                   DYX,DYY,DYZ, &
                                   DZX,DZY,DZZ, &
                                   Elem%Mu, Elem%Lambda, Elem%Cij, &
                                   m1,m2,m3, n_solid, &
                                   Elem%onemSbeta, Elem%sSimu(i_simu)%R_xx_, Elem%sSimu(i_simu)%R_yy_, &
                                   Elem%sSimu(i_simu)%R_xy_, Elem%sSimu(i_simu)%R_xz_, Elem%sSimu(i_simu)%R_yz_)
      call attenuation_update(Elem%sSimu(i_simu)%epsilondev_xx_,Elem%sSimu(i_simu)%epsilondev_yy_, &
                              Elem%sSimu(i_simu)%epsilondev_xy_,Elem%sSimu(i_simu)%epsilondev_xz_,Elem%sSimu(i_simu)%epsilondev_yz_, &
                              epsilondev_xx_loc,epsilondev_yy_loc, &
                              epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                              Elem%sSimu(i_simu)%R_xx_,Elem%sSimu(i_simu)%R_yy_, &
                              Elem%sSimu(i_simu)%R_xy_,Elem%sSimu(i_simu)%R_xz_,Elem%sSimu(i_simu)%R_yz_, &
                              Elem%factor_common_3,Elem%alphaval_3,Elem%betaval_3,Elem%gammaval_3, &
                              Elem%Mu, m1,m2,m3, n_solid)
      deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
   else
      call calcul_forces_aniso(Fox,Foy,Foz,  &
                               Elem%Invgrad(:,:,:,0,0), &
                               Elem%Invgrad(:,:,:,1,0), &
                               Elem%Invgrad(:,:,:,2,0), &
                               Elem%Invgrad(:,:,:,0,1), &
                               Elem%Invgrad(:,:,:,1,1), &
                               Elem%Invgrad(:,:,:,2,1), &
                               Elem%Invgrad(:,:,:,0,2), &
                               Elem%Invgrad(:,:,:,1,2), &
                               Elem%Invgrad(:,:,:,2,2), &
                               htprimex, htprimey, htprimez, &
                               Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                               DXX,DXY,DXZ, &
                               DYX,DYY,DYZ, &
                               DZX,DZY,DZZ, &
                               Elem%Cij, &
                               m1,m2,m3)
   endif
else
   if (n_solid>0) then
      call calcul_forces_att(Fox,Foy,Foz, &
                             Elem%Invgrad(:,:,:,0,0), &
                             Elem%Invgrad(:,:,:,1,0), &
                             Elem%Invgrad(:,:,:,2,0), &
                             Elem%Invgrad(:,:,:,0,1), &
                             Elem%Invgrad(:,:,:,1,1), &
                             Elem%Invgrad(:,:,:,2,1), &
                             Elem%Invgrad(:,:,:,0,2), &
                             Elem%Invgrad(:,:,:,1,2), &
                             Elem%Invgrad(:,:,:,2,2), &
                             htprimex, htprimey, htprimez, &
                             Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                             DXX,DXY,DXZ, &
                             DYX,DYY,DYZ, &
                             DZX,DZY,DZZ, &
                             Elem%Mu, Elem%Lambda, &
                             m1,m2,m3, n_solid, &
                             Elem%onemSbeta, Elem%sSimu(i_simu)%R_xx_, Elem%sSimu(i_simu)%R_yy_, &
                             Elem%sSimu(i_simu)%R_xy_, Elem%sSimu(i_simu)%R_xz_, Elem%sSimu(i_simu)%R_yz_)
      call attenuation_update(Elem%sSimu(i_simu)%epsilondev_xx_,Elem%sSimu(i_simu)%epsilondev_yy_, &
                              Elem%sSimu(i_simu)%epsilondev_xy_,Elem%sSimu(i_simu)%epsilondev_xz_,Elem%sSimu(i_simu)%epsilondev_yz_, &
                              epsilondev_xx_loc,epsilondev_yy_loc, &
                              epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                              Elem%sSimu(i_simu)%R_xx_,Elem%sSimu(i_simu)%R_yy_, &
                              Elem%sSimu(i_simu)%R_xy_,Elem%sSimu(i_simu)%R_xz_,Elem%sSimu(i_simu)%R_yz_, &
                              Elem%factor_common_3,Elem%alphaval_3,Elem%betaval_3,Elem%gammaval_3, &
                              Elem%Mu, m1,m2,m3, n_solid)
      deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
   else
      call calcul_forces(Fox,Foy,Foz,  &
                         Elem%Invgrad(:,:,:,0,0), &
                         Elem%Invgrad(:,:,:,1,0), &
                         Elem%Invgrad(:,:,:,2,0), &
                         Elem%Invgrad(:,:,:,0,1), &
                         Elem%Invgrad(:,:,:,1,1), &
                         Elem%Invgrad(:,:,:,2,1), &
                         Elem%Invgrad(:,:,:,0,2), &
                         Elem%Invgrad(:,:,:,1,2), &
                         Elem%Invgrad(:,:,:,2,2), &
                         htprimex, htprimey, htprimez, &
                         Elem%Jacob, Elem%wgtx, Elem%wgty, Elem%wgtz, &
                         DXX,DXY,DXZ, &
                         DYX,DYY,DYZ, &
                         DZX,DZY,DZZ, &
                         Elem%Mu, Elem%Lambda, &
                         m1,m2,m3)
   endif
endif

Elem%sSimu(i_simu)%Forces(:,:,:,0) = -Fox
Elem%sSimu(i_simu)%Forces(:,:,:,1) = -Foy
Elem%sSimu(i_simu)%Forces(:,:,:,2) = -Foz



return
end subroutine forces_int
