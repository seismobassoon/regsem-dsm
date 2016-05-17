program mesher


use module_ellipticity

implicit none

type :: processor
integer, dimension(:), pointer :: Obj
end type

type :: souvenir
type(processor), dimension(:), pointer :: rank 
end type

type(souvenir), dimension(:), pointer :: memory
doubleprecision :: x_len,y_len,z_len, dx,dy,dz, rayon, R, X, Y, D, ratio, Clat,Clong, alpha, dRmax, &
                   xa,ya,za, xs,ys,zs, theta,phi, limit_inf,limit_med, dist, Ar,Atheta,Aphi, &
                   epsil_xi,epsil_eta,epsil_r, latmin,latmax,dlat, lonmin,lonmax,dlong, Rterre, P2
doubleprecision, parameter :: pi = 3.141592653, seuil = 0.2, vs_surface = 3200.00
doubleprecision, dimension(0:2) :: tmp1D
doubleprecision, dimension(0:3) :: weight, rad_loc
doubleprecision, dimension(0:11) :: rad, vs
doubleprecision, dimension(:), allocatable :: radius, xco,yco,zco, dR, Rsph
doubleprecision, dimension(0:2,0:2) :: rot
doubleprecision, dimension(0:7,0:2) :: coord_vertices
doubleprecision, dimension(:,:), allocatable :: coord_nodes, tmp2D, dzmin, dzmax, &
                                                moho,depth_moho, topo,altitude_topo
integer :: i,j,k, n,m, nel, i_count, num, n_elem, n_elem_xy, n_points, n_pts_xy, nx,ny,nz_tot, &
           edgecut, nparts, proc, ok, n_vertices, n_faces, n_edges, n_nodes, nf,ne,nv, &
           neighbor, neighbor_face, neighbor_edge, nods_per_elem, n_layers, model, nb_couches_model, &
           i1,i2,j1,j2,k1,k2, n_out
integer, parameter :: n_dim = 3, etype = 3, numflag = 0, wgtflag = 2, options = 0
integer, dimension(:), allocatable :: corner, neighbor_corner, counter, tmp1D_int, nz, elmnts, &
                                      dxadj, dxadjncy, part, vwgt, adjwgt, Material, t_reversal, &
                                      moho_position, topo_position, &
                                      nelem_in_proc, which_vertices, nf_shared, ne_shared, nv_shared
integer, dimension(:,:), allocatable :: elmnts_local, which_elem_in_proc, nodes, tmp2D_int, &
                                        faces, faces_shared, edges, edges_shared, vertices_shared
logical :: curve, any_PML, free_surface_pres, topo_log, ellipticity, random_trac, mirror
logical, parameter :: output_chunk = .false.
logical, dimension(:), allocatable :: L_Proc
character(len=1) :: yes_or_not
character(len=13) :: meshfilename


!!! Inputs !!!

write (*,*) "Do you want to take into account the Earth sphericity ? [y/n]"
read (*,*) yes_or_not
if (yes_or_not == "y") then
   curve = .true.
   write (*,*) "Introduce the coordinates (lat,long) of the center"
   read (*,*) Clat, Clong
   write (*,*) "Introduce the lengths of the horizontal edges: Xi and Eta (max: 90 degrees)"
   read (*,*) x_len, y_len
   write (*,*) "Introduce the horizontal space step:"
   read (*,*) dx;   dy = dx
   write (*,*) "Do you want PREM interfaces (1)"
   write (*,*) "            PREM interfaces + 3D Moho (2)" ! Les 3 premieres interfaces de PREM sont alors virees
   write (*,*) "            PREM interfaces + 3D Moho - 220km interface (3)" ! Les 4 premieres interfaces de PREM sont alors virees
   write (*,*) "            another 1D model (4)"
   write (*,*) "            SVEMum_SAW24B16LM (5)"
   write (*,*) "            PREM interfaces - 15km interface (6)"
   read (*,*) model
   if (model==4) then
      ! Ici on peut construire une boule a couches de rayon quelconque
      ! Une topo de surface pourra etre prise en compte
      ! Les oceans pourront aussi etre consideres. Il faudra alors ajuster le parametre Rterre dans angles.f90
      ! Les vitesses du modele qui sera associe a ce maillage devront croitre avec la profondeur
      write (*,*) "Introduce the number of layers"
      read (*,*) n_layers
      write (*,"(a11,i2.2,a18)") " Introduce ", n_layers+1, " radii (in meter)"
      allocate (radius(0:n_layers))
      read (*,*) (radius(i),i=0,n_layers)
      Rterre = radius(n_layers)
      ellipticity = .false.
   else
      if ((model==2 .or. model==3) .and. (dx<0.199999))then
          print *,"The horizontal step has to be larger than 0.2 degree."
          print *,"Otherwise you'll probably have an under-sampled crust."
          stop ! Sinon on aurait plusieurs couches d'elements dans la croute
      endif
      call define_model(model,nb_couches_model,rad,vs)
      Rterre = 6371000.d0
      rad(nb_couches_model-1) = Rterre
      vs(nb_couches_model-1) = vs_surface
      write (*,*) "Introduce the depth (in meter) of the bottom"
      read (*,*) rayon
      rayon = Rterre - rayon
      bottom : do i = 0,nb_couches_model-1
         if (rayon<rad(i)) then
            n_layers = nb_couches_model - i
            exit bottom
         endif
      enddo bottom
      call bottom_remarks(rayon,rad,i,dx,Rterre,model)
      allocate (radius(0:n_layers))
      radius(0) = rayon
      do j = 1,n_layers
         radius(j) = rad(i)
         i = i + 1
      enddo
      write (*,*) "Do you want to take into account the Earth ellipticity ? [y/n]"
      read (*,*) yes_or_not
      ellipticity = .false.
      if (yes_or_not == "y") then
          ellipticity = .true.
          call init_ellipticity()
      endif
   endif
   write (*,*) "Do you want a topography at the surface ? [y/n]"
   read (*,*) yes_or_not
   topo_log = .false.
   if (yes_or_not == "y") topo_log = .true.
   write (*,*) "How many nodes do you want for the geometry of each element: 8, 20 or 27 ?"
   read (*,*) nods_per_elem
else
   curve = .false.
   nods_per_elem = 8
   model = 0
   ellipticity = .false.
   topo_log = .false.
   write (*,*) "Introduce the lengths of the edges: x_len, y_len, z_len (in meter)"
   read (*,*) x_len, y_len, z_len
   write (*,*) "Introduce the space steps: dx, dy, dz"
   read (*,*) dx, dy, dz
endif

write (*,*) "Do you want PML ? [y/n]"
read (*,*) yes_or_not
any_PML = .false.
if (yes_or_not == "y") then
   any_PML = .true.
   write (*,*) "Introduce the polynomial order for a PML element, and then for a normal element"
   read (*,*) i, j;   ratio = (i+1.)/(j+1.)
   write (*,*) "Do you want a free surface at the top ? [y/n]"
   read (*,*) yes_or_not
   free_surface_pres = .false.
   if (yes_or_not == "y") free_surface_pres = .true.
endif

write (*,*) "Introduce the number of parts to partition the mesh"
read (*,*) nparts

random_trac = .false.
if (curve .and. nods_per_elem==27) then
   write (*,*) "Do you want a file for random tractions ? [y/n]"
   read (*,*) yes_or_not
   if (yes_or_not == "y") then
      random_trac = .true.
      write (*,*) "Introduce the direction of the tractions : A_r, A_theta, A_phi"
      read (*,*) Ar, Atheta, Aphi
      open (25,file="random_trac.table")
      write (25,*) Ar, Atheta, Aphi
   endif
endif

mirror = .false.
if (any_PML) then
   write (*,*) "Do you want a time-reversal mirror ? [y/n]"
   read (*,*) yes_or_not
   if (yes_or_not == "y") then
      mirror = .true.
      write (*,*) "Define the mirror with the number of elements that you"
      if (free_surface_pres) then
         write (*,*) "want to remove in each direction (-x,x,-y,y,-z) ?"
         read (*,*) i1,i2,j1,j2,k1;
         i1 = i1+1;   i2 = i2+1;   j1 = j1+1;   j2 = j2+1;   k1 = k1+1;   k2 = 0
      else
         write (*,*) "want to remove in each direction (-x,x,-y,y,-z,z) ?"
         read (*,*) i1,i2,j1,j2,k1,k2
         i1 = i1+1;   i2 = i2+1;   j1 = j1+1;   j2 = j2+1;   k1 = k1+1;   k2 = k2+1
      endif
   endif
endif


!!! Defining the number of elements !!!

nx = int (x_len/dx)
ny = int (y_len/dy)

if (nx*dx < x_len-dx/2.) nx = nx+1
if (ny*dy < y_len-dy/2.) ny = ny+1

dx = x_len/nx
dy = y_len/ny

if (curve) then
   dRmax = pi*dx*radius(n_layers)/180.d0
   allocate (nz(0:n_layers))
   allocate (dR(0:n_layers))
   nz_tot = 0
   do i = 0,n_layers-1
      z_len = radius(i+1) - radius(i)
      if (model==4) then
         dz = dRmax
      else
         dz = dRmax * vs(nb_couches_model-(n_layers-i))/vs(nb_couches_model-1)
         if (egal(dz,0.d0,1.d0)==.true.) dz = dRmax ! Quand on se trouve dans le noyau
      endif
      X = z_len/dz
      Y = X - aint(X)
      if (Y/X > seuil) then
         nz(i) = int(X) + 1
      else
         nz(i) = int(X)
      endif
      dR(i) = z_len/nz(i)
      nz_tot = nz_tot + nz(i)
   enddo
   nz(n_layers) = 1
   dR(n_layers) = 0
else
   nz_tot = int (z_len/dz)
   if (nz_tot*dz < z_len-dz/2.) nz_tot = nz_tot+1
   dz = z_len/nz_tot
endif

n_points = (nx+1) * (ny+1) * (nz_tot+1)
n_pts_xy = (nx+1) * (ny+1)
n_elem = nx * ny * nz_tot
n_elem_xy = nx * ny


!!! Associating to each element 8 vertices by using a global numbering !!!
allocate (elmnts(0:8*n_elem-1))
do nel = 0, n_elem-1 
    k = nel / n_elem_xy
    j = nel - k*n_elem_xy; j = j/nx
    i = nel - k*n_elem_xy - j*nx 
    i_count = i + j*(nx+1) + k*n_pts_xy
    elmnts(nel*8:nel*8+7) = (/ i_count, i_count+1, i_count+nx+2, i_count+nx+1, i_count+n_pts_xy, &
    i_count+n_pts_xy+1, i_count+n_pts_xy+nx+2, i_count+n_pts_xy+nx+1 /)
enddo

!!! Defining the coordinates of each vertex !!!
allocate (xco(0:n_points-1))
allocate (yco(0:n_points-1))
allocate (zco(0:n_points-1))
if (curve) then
   y_len = -y_len/2.d0
   x_len = -x_len/2.d0
   i_count = 0
   do n = 0,n_layers
      rayon = radius(n)
      dz = dR(n)
      do k = 0,nz(n)-1
         R = rayon + k*dz
         do j = 0,ny
            Y = y_len + j*dy
            do i = 0,nx
               X = x_len + i*dx
               xco(i_count) = X
               yco(i_count) = Y
               zco(i_count) = R
               i_count = i_count + 1
            enddo
         enddo
      enddo
   enddo
else
   do k = 0,nz_tot
      do j = 0,ny
         do i = 0,nx
            i_count = i + j*(nx+1) + k*n_pts_xy
            xco(i_count) = i * dx
            yco(i_count) = j * dy
            zco(i_count) = k * dz
         enddo
      enddo
   enddo
endif

!!! Defining the rotation matrix !!!
if (curve) then
   Clat = 90.d0 - Clat
   if (Clong<0) Clong = 360.d0 + Clong
   alpha = 0
   Clat = deg2rad(Clat); Clong = deg2rad(Clong); alpha = deg2rad(alpha);
   rot(0,0) = dcos(alpha)*dcos(Clat)*dcos(Clong)-dsin(alpha)*dsin(Clong)
   rot(0,1) = -dsin(alpha)*dcos(Clat)*dcos(Clong)-dcos(alpha)*dsin(Clong)
   rot(0,2) = dsin(Clat)*dcos(Clong)
   rot(1,0) = dcos(alpha)*dcos(Clat)*dsin(Clong)+dsin(alpha)*dcos(Clong)
   rot(1,1) = -dsin(alpha)*dcos(Clat)*dsin(Clong)+dcos(alpha)*dcos(Clong)
   rot(1,2) = dsin(Clat)*dsin(Clong)
   rot(2,0) = -dcos(alpha)*dsin(Clat)
   rot(2,1) = dsin(alpha)*dsin(Clat)
   rot(2,2) = dcos(Clat)
endif

!!! Defining the topography of the Moho !!!
!!! On agit comme s'il y avait 27 ctrl_pts car c'est le cas general (8 et 20 sont des sous-cas de 27) !!!
if (model==2 .or. model==3) then
    open (11, file="Moho.asc", status="old")
    read (11,*) latmin, latmax, dlat
    if (mod(latmax-latmin,dlat)/=0) then
        print *,"In the file Moho.asc the latitude step doesn't fit the latitude limits"
        stop
    endif
    read (11,*) lonmin, lonmax, dlong
    if (lonmin<0.d0)   lonmin = 360.d0 + lonmin
    if (lonmax<0.d0)   lonmax = 360.d0 + lonmax
    if (mod(lonmax-lonmin,dlong)/=0) then
        print *,"In the file Moho.asc the longitude step doesn't fit the longitude limits"
        stop
    endif
    m = int((latmax-latmin)/dlat) + 1
    n = int((lonmax-lonmin)/dlong) + 1
    allocate (moho(0:m-1,0:n-1))
    do i = 0,m-1
     do j = 0,n-1
         read (11,*) moho(i,j)
     enddo
    enddo
    moho(:,:) = moho(:,:)*1000.d0
    allocate (depth_moho(0:2*nx,0:2*ny))
    open (17,file="moho.out")
    do j = 0,2*ny
        Y = y_len + j*dy/2.d0
        do i = 0,2*nx
            X = x_len + i*dx/2.d0
            ! Passage en coordonnees cartesiennes
            xa = tan (pi*X/180.d0)
            ya = tan (pi*Y/180.d0)
            D = sqrt(1.d0 + xa**2 + ya**2)
            xa = xa/D;   ya = ya/D;   za = 1/D
            ! Faire la rotation du chunk de ref vers le chunk reel
            xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
            ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
            zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
            ! Passage en spherique
            call cart2sph(xs,ys,zs,alpha,theta,phi)
            ! Profondeur du Moho en (theta,phi)
            call read_moho(theta,phi,moho,m,n,latmin,latmax,dlat,lonmin,lonmax,dlong,depth_moho(i,j))
            write(17,*) 90.d0-theta*180.d0/pi,phi*180.d0/pi,depth_moho(i,j)
        enddo
    enddo
    close(17)
    deallocate (moho)
    allocate (dzmin(0:1,0:nparts-1));   dzmin = 1000000
    allocate (dzmax(0:1,0:nparts-1));   dzmax = 0
    limit_inf = radius(n_layers-2) + dR(n_layers-2)/10.d0
    limit_med = radius(n_layers-1) - dR(n_layers-2)/10.d0
endif

!!! Defining the topography of the surface!!!
!!! On procede comme pour le Moho !!!
if (topo_log) then
    open (11, file="Topo.asc", status="old")
    read (11,*) latmin, latmax, dlat
    if (mod(latmax-latmin,dlat)/=0) then
        print *,"In the file Topo.asc the latitude step doesn't fit the latitude limits"
        stop
    endif
    read (11,*) lonmin, lonmax, dlong
    if (lonmin<0.d0)   lonmin = 360.d0 - lonmin
    if (lonmax<0.d0)   lonmax = 360.d0 - lonmax
    if (mod(lonmax-lonmin,dlong)/=0) then
        print *,"In the file Topo.asc the longitude step doesn't fit the longitude limits"
        stop
    endif
    m = int((latmax-latmin)/dlat) + 1
    n = int((lonmax-lonmin)/dlong) + 1
    allocate (topo(0:m-1,0:n-1))
    do i = 0,m-1
     do j = 0,n-1
         read (11,*) topo(i,j)
     enddo
    enddo
    allocate (altitude_topo(0:2*nx,0:2*ny))
    open (17,file="topo.out")
    do j = 0,2*ny
        Y = y_len + j*dy/2.d0
        do i = 0,2*nx
            X = x_len + i*dx/2.d0
            ! Passage en coordonnees cartesiennes
            xa = tan (pi*X/180.d0)
            ya = tan (pi*Y/180.d0)
            D = sqrt(1.d0 + xa**2 + ya**2)
            xa = xa/D;   ya = ya/D;   za = 1/D
            ! Faire la rotation du chunk de ref vers le chunk reel
            xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
            ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
            zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
            ! Passage en spherique
            call cart2sph(xs,ys,zs,alpha,theta,phi)
            ! Topo en (theta,phi). On utilise la routine ecrite pour le Moho.
            call read_moho(theta,phi,topo,m,n,latmin,latmax,dlat,lonmin,lonmax,dlong,altitude_topo(i,j))
            write(17,*) 90.d0-theta*180.d0/pi,phi*180.d0/pi,altitude_topo(i,j)
        enddo
    enddo
    close(17)
    deallocate (topo)
    limit_med = radius(n_layers-1) - dR(n_layers-2)/10.d0
endif

!!! Signing the elements above and below the Moho !!!
!!! On suppose ici qu'il n'y a qu'une seule couche d'elements entre le Moho et la surface !!!
allocate (moho_position(0:n_elem-1))
moho_position = 0
if (model==2 .or. model==3) then
    do k = 0,nz_tot-1
      do j = 0,ny-1
        do i = 0,nx-1
            i_count = i + j*nx + k*n_elem_xy
            if (k==nz_tot-2) then
                moho_position(i_count) = -1
            else if (k==nz_tot-1) then
                moho_position(i_count) = 1
            endif
        enddo
      enddo
    enddo
endif

!!! Signing the elements below the 3D surface !!!
allocate (topo_position(0:n_elem-1))
topo_position = 0
if (topo_log) then
    do k = 0,nz_tot-1
      do j = 0,ny-1
        do i = 0,nx-1
            i_count = i + j*nx + k*n_elem_xy
            if (k==nz_tot-1) topo_position(i_count) = 1
        enddo
      enddo
    enddo
endif

!!! Signing the elements that define the time-reversal mirror !!!
allocate (t_reversal(0:n_elem-1))
t_reversal = 0
if (mirror) then
    k2 = nz_tot-1-k2;   j2 = ny-1-j2;   i2 = nx-1-i2
    do k = k1,k2
      do j = j1,j2
        do i = i1,i2
            i_count = i + j*nx + k*n_elem_xy

            if (k==k1) then
                t_reversal(i_count) = 1
            else if (k==k2 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 6
            else if (j==j1) then
                t_reversal(i_count) = 2
            else if (j==j2) then
                t_reversal(i_count) = 4
            else if (i==i1) then
                t_reversal(i_count) = 5
            else if (i==i2) then
                t_reversal(i_count) = 3
            endif

            if (k==k1 .and. j==j1) then
                t_reversal(i_count) = 7
            else if (k==k1 .and.i==i2) then
                t_reversal(i_count) = 8
            else if (k==k1 .and. j==j2) then
                t_reversal(i_count) = 9
            else if (k==k1 .and. i==i1) then
                t_reversal(i_count) = 10
            else if (i==i1 .and. j==j1) then
                t_reversal(i_count) = 11
            else if (i==i2 .and. j==j1) then
                t_reversal(i_count) = 12
            else if (i==i2 .and. j==j2) then
                t_reversal(i_count) = 13
            else if (i==i1 .and. j==j2) then
                t_reversal(i_count) = 14
            else if (k==k2 .and. j==j1 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 15
            else if (k==k2 .and. i==i2 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 16
            else if (k==k2 .and. j==j2 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 17
            else if (k==k2 .and. i==i1 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 18
            endif

            if (k==k1 .and. j==j1 .and. i ==0) then
                t_reversal(i_count) = 19
            else if (k==k1 .and. i==i2 .and. j==j1) then
                t_reversal(i_count) = 20
            else if (k==k1 .and. j==j2 .and. i==i2) then
                t_reversal(i_count) = 21
            else if (k==k1 .and. i==i1 .and. j==j2) then
                t_reversal(i_count) = 22
            else if (k==k2 .and. j==j1 .and. i==i1 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 23
            else if (k==k2 .and. i==i2 .and. j==j1 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 24
            else if (k==k2 .and. j==j2 .and. i==i2 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 25
            else if (k==k2 .and. i==i1 .and. j==j2 .and. (.not.free_surface_pres)) then
                t_reversal(i_count) = 26
            endif

        enddo
      enddo
    enddo
endif

!!! To each element we associate a number which refers to a material !!!
allocate (Material(0:n_elem-1))
Material = 0
if (any_PML) then
    do k = 0,nz_tot-1
      do j = 0,ny-1
        do i = 0,nx-1
            i_count = i + j*nx + k*n_elem_xy

            if (k==0) then
                Material(i_count) = 1
            else if (k==nz_tot-1 .and. (.not.free_surface_pres)) then
                Material(i_count) = 6
            else if (j==0) then
                Material(i_count) = 2
            else if (j==ny-1) then
                Material(i_count) = 4
            else if (i==0) then
                Material(i_count) = 5
            else if (i==nx-1) then
                Material(i_count) = 3
            endif

            if (k==0 .and. j==0) then
                Material(i_count) = 7
            else if (k==0 .and.i==nx-1) then
                Material(i_count) = 8
            else if (k==0 .and. j==ny-1) then
                Material(i_count) = 9
            else if (k==0 .and. i==0) then
                Material(i_count) = 10
            else if (i==0 .and. j==0) then
                Material(i_count) = 11
            else if (i==nx-1 .and. j==0) then
                Material(i_count) = 12
            else if (i==nx-1 .and. j==ny-1) then
                Material(i_count) = 13
            else if (i==0 .and. j==ny-1) then
                Material(i_count) = 14
            else if (k==nz_tot-1 .and. j==0 .and. (.not.free_surface_pres)) then
                Material(i_count) = 15
            else if (k==nz_tot-1 .and. i==nx-1 .and. (.not.free_surface_pres)) then
                Material(i_count) = 16
            else if (k==nz_tot-1 .and. j==ny-1 .and. (.not.free_surface_pres)) then
                Material(i_count) = 17
            else if (k==nz_tot-1 .and. i==0 .and. (.not.free_surface_pres)) then
                Material(i_count) = 18
            endif

            if (k==0 .and. j==0 .and. i==0) then
                Material(i_count) = 19
            else if (k==0 .and. i==nx-1 .and. j==0) then
                Material(i_count) = 20
            else if (k==0 .and. j==ny-1 .and. i==nx-1) then
                Material(i_count) = 21
            else if (k==0 .and. i==0 .and. j==ny-1) then
                Material(i_count) = 22
            else if (k==nz_tot-1 .and. j==0 .and. i==0 .and. (.not.free_surface_pres)) then
                Material(i_count) = 23
            else if (k==nz_tot-1 .and. i==nx-1 .and. j==0 .and. (.not.free_surface_pres)) then
                Material(i_count) = 24
            else if (k==nz_tot-1 .and. j==ny-1 .and. i==nx-1 .and. (.not.free_surface_pres)) then
                Material(i_count) = 25
            else if (k==nz_tot-1 .and. i==0 .and. j==ny-1 .and. (.not.free_surface_pres)) then
                Material(i_count) = 26
            endif

        enddo
      enddo
    enddo
endif
allocate (vwgt(0:n_elem-1))
do i = 0,3
    weight(i) = anint(ratio**i)
enddo
do nel = 0, n_elem-1
    if (Material(nel)<1) then
        vwgt(nel) = weight(0)
    else if (Material(nel)<7) then
        vwgt(nel) = weight(1)
    else if (Material(nel)<19) then
        vwgt(nel) = weight(2)
    else
        vwgt(nel) = weight(3)
    endif
enddo

!!! Partitioning the domain by using METIS !!!
allocate (dxadj(0:n_elem))
allocate (dxadjncy(0:8*n_elem-1))
call METIS_MeshToDual (n_elem, n_points, elmnts, etype, numflag, dxadj, dxadjncy)
allocate (adjwgt(0:dxadj(n_elem)-1)); adjwgt = 0
allocate (part(0:n_elem-1))
if (nparts==1) then
    do nel = 0,n_elem-1
        part(nel) = 0
    enddo
else
    call METIS_PartGraphKway (n_elem, dxadj, dxadjncy(0:dxadj(n_elem)-1), vwgt, adjwgt, &
                              wgtflag, numflag, nparts, options, edgecut, part)
endif
!!! Now there will be two different numberings: !!!
!!! a local one (i.e. relative to a processor) and a global one (i.e. which concerns the whole domain) !!!

!!! Defining the number of elements per processor !!!
allocate (nelem_in_proc(0:nparts-1))
nelem_in_proc = 0
do nel = 0,n_elem-1
    nelem_in_proc(part(nel)) = nelem_in_proc(part(nel)) + 1
enddo

!!! Defining for each processor which elements are inside !!!
!!! "Counter" refers to the local numberings and "nel" to the global numbering !!!
!!! Note that the elements of a processor are naturally sorted according to the global numbering !!!
allocate (counter(0:nparts-1))
counter = 0
allocate (which_elem_in_proc(0:nparts-1,0:maxval(nelem_in_proc)-1))
do nel = 0,n_elem-1
    num = part(nel)
    which_elem_in_proc(num,counter(num)) = nel
    counter(num) = counter(num) + 1
enddo
deallocate (counter)

!!! Allocating "memory" !!!
!!! Allow to establish the correspondence between objects shared by different processors !!!
allocate (memory(0:n_elem-1))
do nel = 0,n_elem-1
    if (part(nel) /= nparts-1) then
        allocate (memory(nel)%rank(part(nel)+1:nparts-1))
        do proc = part(nel)+1,nparts-1
            allocate (memory(nel)%rank(proc)%Obj(0:17))
        enddo
    endif
enddo


!!! LET'S CONSIDER A PROCESSOR AND CREATE THE FIELDS REQUIRED TO BUILD ITS MESHFILE !!!
call system ("rm -f mesh4spec.???")
meshfilename(1:10) = "mesh4spec."
do proc = 0,nparts-1

 !!! Defining the vertices which belong to the processor !!!
 !!! These vertices will be sorted according to the global numbering !!!
 allocate (which_vertices(0:8*nelem_in_proc(proc)-1))
 n_vertices = 0
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     do i = 0,7
         num = elmnts(8*nel+i)
         ok = 1
         check : do j = 0,n_vertices-1
             if (which_vertices(j)==num) then
                 ok = 0
                 exit check
             endif
         enddo check
         if (ok==1) then
             which_vertices(n_vertices) = num
             n_vertices = n_vertices+1
         endif
     enddo
 enddo
 allocate (tmp1D_int(0:n_vertices-1))
 tmp1D_int(0:n_vertices-1) = which_vertices(0:n_vertices-1)
 deallocate (which_vertices)
 allocate (which_vertices(0:n_vertices-1))
 which_vertices(0:n_vertices-1) = tmp1D_int(0:n_vertices-1)
 deallocate (tmp1D_int)
 call sort(which_vertices,n_vertices)

 !!! Associating to each element 8 vertices by using a local numbering !!!
 allocate (elmnts_local(0:nelem_in_proc(proc)-1,0:7))
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     do i = 0,7
         num = elmnts(8*nel+i)
         build : do j = 0,n_vertices-1
             if (which_vertices(j)==num) then
                 elmnts_local(n,i) = j
                 exit build
             endif
         enddo build
     enddo
 enddo

 !!! Associating to each element 6 faces by using a local numbering !!!
 !!! Defining the faces shared with another processor !!!
 allocate (faces(0:nelem_in_proc(proc)-1,0:5))
 allocate (faces_shared(0:nparts-1,0:6*nelem_in_proc(proc)-1))
 allocate (nf_shared(0:nparts-1))
 nf_shared = 0
 allocate (corner(0:3))
 allocate (neighbor_corner(0:3))
 n_faces = 0
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering
     do nf = 0,5   ! nf indicates which face of the element we're considering
         select case (nf)   ! Here we pick up the vertices (in the global numbering) which define the face
         case (0)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          corner(2) = elmnts(8*nel+2)
          corner(3) = elmnts(8*nel+3)
         case (1)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
          corner(2) = elmnts(8*nel+4)
          corner(3) = elmnts(8*nel+5)
         case (2)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+2)
          corner(2) = elmnts(8*nel+5)
          corner(3) = elmnts(8*nel+6)
         case (3)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+3)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+7)
         case (4)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+3)
          corner(2) = elmnts(8*nel+4)
          corner(3) = elmnts(8*nel+7)
         case (5)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+5)
          corner(2) = elmnts(8*nel+6)
          corner(3) = elmnts(8*nel+7)
         end select
         find0 : do i = dxadj(nel), dxadj(nel+1)-1   ! We look at a neighbor of the element
             neighbor = dxadjncy(i)
             find1 : do j = 0,3   ! DOES THE FACE BELONG TO THIS NEIGHBOR ?
                 num = corner(j)
                 ok = 0
                 find2 : do k = 0,7
                     if (elmnts(8*neighbor+k)==num) then
                         neighbor_corner(j) = k
                         ok = 1
                         exit find2
                     endif
                 enddo find2
                 if (ok==0) exit find1   ! NO, so let's see another neighbor
                 if (j==3) then   ! YES
                     call sort(neighbor_corner,4)
                     if (neighbor_corner(0)==0) then   ! So which face of the neighbor is it ?
                     if (neighbor_corner(3)==3) neighbor_face = 0
                     if (neighbor_corner(3)==5) neighbor_face = 1
                     if (neighbor_corner(3)==7) neighbor_face = 4
                     else if (neighbor_corner(0)==1) then
                         neighbor_face = 2
                     else if (neighbor_corner(0)==2) then
                         neighbor_face = 3
                     else if (neighbor_corner(0)==4) then
                         neighbor_face = 5
                     else
                         print *,"Coherency Pb between faces and vertices of an element"
                     endif
                     if (part(neighbor)==proc) then   ! The neighbor is on the same processor than the element
                         if (neighbor>nel) then   ! The neighbor is an element we've never seen
                             faces(n,nf) = n_faces
                             n_faces = n_faces + 1
                         else   ! The neighbor is an element we've ever seen
                             g2l : do i_count = 0,n-1
                                 if (which_elem_in_proc(proc,i_count)==neighbor) then
                                     faces(n,nf) = faces(i_count,neighbor_face)
                                     exit g2l
                                 endif
                             enddo g2l
                         endif   
                     else   ! The neighbor is not on the same processor than the element
                         faces(n,nf) = n_faces
                         num = part(neighbor)
                         !!! Ensuring the correspondence between the faces shared !!!
                         if (num<proc) then   ! We've ever seen the processor of the neighbor
                             faces_shared(num,memory(neighbor)%rank(proc)%Obj(neighbor_face)) = n_faces
                         else   ! We've never seen the processor of the neighbor
                             faces_shared(num,nf_shared(num)) = n_faces
                             memory(nel)%rank(num)%Obj(nf) = nf_shared(num)
                         endif
                         nf_shared(num) = nf_shared(num) + 1
                         n_faces = n_faces + 1
                     endif
                     exit find0
                 endif
             enddo find1
         enddo find0
         if (ok==0) then   ! The face is not shared by a neighbor
             faces(n,nf) = n_faces
             n_faces = n_faces + 1
         endif
     enddo
 enddo
 deallocate (corner, neighbor_corner)
 allocate (tmp2D_int(0:nparts-1,0:maxval(nf_shared)-1))
 tmp2D_int(0:nparts-1,0:maxval(nf_shared)-1) = faces_shared(0:nparts-1,0:maxval(nf_shared)-1)
 deallocate (faces_shared)
 allocate (faces_shared(0:nparts-1,0:maxval(nf_shared)-1))
 faces_shared(0:nparts-1,0:maxval(nf_shared)-1) = tmp2D_int(0:nparts-1,0:maxval(nf_shared)-1)
 deallocate (tmp2D_int)

 !!! Associating to each element 12 edges by using a local numbering !!!
 !!! Defining the edges shared with others processors !!!
 allocate (edges(0:nelem_in_proc(proc)-1,0:11))
 allocate (L_Proc(0:nparts-1))
 allocate (edges_shared(0:nparts-1,0:12*nelem_in_proc(proc)-1))
 allocate (ne_shared(0:nparts-1))
 ne_shared = 0
 allocate (corner(0:1))
 allocate (neighbor_corner(0:1))
 n_edges = 0
 do n = 0,nelem_in_proc(proc)-1   ! n is the number of the considered element in the local numbering
     nel = which_elem_in_proc(proc,n)   ! nel is the number of the considered element in the global numbering
     do ne = 0,11   ! ne indicates which edge of the element we're considering
         select case (ne)   ! Here we pick up the vertices (in the global numbering) which define the edge
         case (0)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+1)
         case (1)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+2)
         case (2)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+3)
         case (3)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+3)
         case (4)
          corner(0) = elmnts(8*nel+1)
          corner(1) = elmnts(8*nel+5)
         case (5)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+5)
         case (6)
          corner(0) = elmnts(8*nel)
          corner(1) = elmnts(8*nel+4)
         case (7)
          corner(0) = elmnts(8*nel+2)
          corner(1) = elmnts(8*nel+6)
         case (8)
          corner(0) = elmnts(8*nel+5)
          corner(1) = elmnts(8*nel+6)
         case (9)
          corner(0) = elmnts(8*nel+6)
          corner(1) = elmnts(8*nel+7)
         case (10)
          corner(0) = elmnts(8*nel+3)
          corner(1) = elmnts(8*nel+7)
         case (11)
          corner(0) = elmnts(8*nel+4)
          corner(1) = elmnts(8*nel+7)
         end select
         findbis0 : do i = 0,n-1   ! HAVE WE EVER SEEN THE EDGE IN ANOTHER ELEMENT OF THE SAME PROCESSOR BEFORE?
             neighbor = which_elem_in_proc(proc,i)
             findbis1 : do j = 0,1
                 num = corner(j)
                 ok = 0
                 findbis2 : do k = 0,7
                     if (elmnts(8*neighbor+k)==num) then
                         neighbor_corner(j) = k
                         ok = 1
                         exit findbis2
                     endif
                 enddo findbis2
                 if (ok==0) exit findbis1
                 if (j==1) then   ! YES
                     call sort(neighbor_corner,2)
                     if (neighbor_corner(0)==0) then
                      if (neighbor_corner(1)==1) neighbor_edge = 0
                      if (neighbor_corner(1)==3) neighbor_edge = 3
                      if (neighbor_corner(1)==4) neighbor_edge = 6
                     else if (neighbor_corner(0)==1) then
                      if (neighbor_corner(1)==2) neighbor_edge = 1
                      if (neighbor_corner(1)==5) neighbor_edge = 4
                     else if (neighbor_corner(0)==2) then
                      if (neighbor_corner(1)==3) neighbor_edge = 2
                      if (neighbor_corner(1)==6) neighbor_edge = 7
                     else if (neighbor_corner(0)==3) then
                      neighbor_edge = 10
                     else if (neighbor_corner(0)==4) then
                      if (neighbor_corner(1)==5) neighbor_edge = 5
                      if (neighbor_corner(1)==7) neighbor_edge = 11
                     else if (neighbor_corner(0)==5) then
                      neighbor_edge = 8
                     else if (neighbor_corner(0)==6) then
                      neighbor_edge = 9
                     else
                      print *,"Coherency Pb between edges and vertices of an element"
                     endif
                     edges(n,ne) = edges(i,neighbor_edge)
                     exit findbis0
                 endif
             enddo findbis1
         enddo findbis0
         if (ok==0 .or. n==0) then   ! NO
             edges(n,ne) = n_edges
             n_edges = n_edges + 1
             L_Proc = .true.
             L_Proc(proc) = .false.
             do i = 0,n_elem-1   ! IS THE EDGE SHARED WITH AN ELEMENT FROM ANOTHER PROCESSOR?
                 if (L_Proc(part(i))) then
                     findter1 : do j = 0,1
                         num = corner(j)
                         ok = 0
                         findter2 : do k = 0,7
                             if (elmnts(8*i+k)==num) then
                                 neighbor_corner(j) = k
                                 ok = 1
                                 exit findter2
                             endif
                         enddo findter2
                         if (ok==0) exit findter1
                         if (j==1) then   ! YES
                             num = part(i)
                             !!! Ensuring the correspondence between the edges shared !!!
                             if (num<proc) then   ! It deals with a processor we've ever seen
                                 call sort(neighbor_corner,2)
                                 if (neighbor_corner(0)==0) then
                                  if (neighbor_corner(1)==1) neighbor_edge = 0
                                  if (neighbor_corner(1)==3) neighbor_edge = 3
                                  if (neighbor_corner(1)==4) neighbor_edge = 6
                                 else if (neighbor_corner(0)==1) then
                                  if (neighbor_corner(1)==2) neighbor_edge = 1
                                  if (neighbor_corner(1)==5) neighbor_edge = 4
                                 else if (neighbor_corner(0)==2) then
                                  if (neighbor_corner(1)==3) neighbor_edge = 2
                                  if (neighbor_corner(1)==6) neighbor_edge = 7
                                 else if (neighbor_corner(0)==3) then
                                  neighbor_edge = 10
                                 else if (neighbor_corner(0)==4) then
                                  if (neighbor_corner(1)==5) neighbor_edge = 5
                                  if (neighbor_corner(1)==7) neighbor_edge = 11
                                 else if (neighbor_corner(0)==5) then
                                  neighbor_edge = 8
                                 else if (neighbor_corner(0)==6) then
                                  neighbor_edge = 9
                                 else
                                  print *,"Coherency Pb between edges and vertices of an element"
                                 endif
                                 edges_shared(num,memory(i)%rank(proc)%Obj(neighbor_edge+6)) = edges(n,ne)
                             else   ! It deals with a processor we've never seen
                                 edges_shared(num,ne_shared(num)) = edges(n,ne)
                                 memory(nel)%rank(num)%Obj(ne+6) = ne_shared(num)
                             endif
                             ne_shared(num) = ne_shared(num) + 1
                             L_Proc(num) = .false.
                         endif
                     enddo findter1
                 endif
             enddo
         endif
     enddo
 enddo
 deallocate (corner, neighbor_corner)
 allocate (tmp2D_int(0:nparts-1,0:maxval(ne_shared)-1))
 tmp2D_int(0:nparts-1,0:maxval(ne_shared)-1) = edges_shared(0:nparts-1,0:maxval(ne_shared)-1)
 deallocate (edges_shared)
 allocate (edges_shared(0:nparts-1,0:maxval(ne_shared)-1))
 edges_shared(0:nparts-1,0:maxval(ne_shared)-1) = tmp2D_int(0:nparts-1,0:maxval(ne_shared)-1)
 deallocate (tmp2D_int)

 !!! Defining vertices shared with others processors !!!
 allocate (vertices_shared(0:nparts-1,0:7*nelem_in_proc(proc)-1))
 allocate (nv_shared(0:nparts-1))
 nv_shared = 0
 do n = 0,n_vertices-1
     L_Proc = .true.
     L_Proc(proc) = .false.
     do i = 0,n_elem-1
         if (L_Proc(part(i))) then
             do k = 0,7
                 if (elmnts(8*i+k)==which_vertices(n)) then
                     num = part(i)
                     !!! The local numbering of the vertices follows the global numbering !!!
                     !!! Using memory to ensure the correspondence between the processors is therefore useless !!!
                     vertices_shared(num,nv_shared(num)) = n
                     nv_shared(num) = nv_shared(num) + 1
                     L_Proc(num) = .false.
                 endif
             enddo
         endif
     enddo
 enddo
 deallocate (L_Proc)
 allocate (tmp2D_int(0:nparts-1,0:maxval(nv_shared)-1))
 tmp2D_int(0:nparts-1,0:maxval(nv_shared)-1) = vertices_shared(0:nparts-1,0:maxval(nv_shared)-1)
 deallocate (vertices_shared)
 allocate (vertices_shared(0:nparts-1,0:maxval(nv_shared)-1))
 vertices_shared(0:nparts-1,0:maxval(nv_shared)-1) = tmp2D_int(0:nparts-1,0:maxval(nv_shared)-1)
 deallocate (tmp2D_int)

 !!! Defining the nodes in the processor !!!
 allocate (coord_nodes(0:nods_per_elem*nelem_in_proc(proc)-1,0:2))
 n_nodes = 0
 do n = 0,nelem_in_proc(proc)-1
     do i = 0,7
         nv = elmnts_local(n,i)
         num = which_vertices(nv)
         coord_vertices(i,0) = xco(num)
         coord_vertices(i,1) = yco(num)
         coord_vertices(i,2) = zco(num)
         num = n_nodes
         call check_vertices(nv,n,elmnts_local(0:n-1,0:7),coord_nodes(num,0:2),coord_vertices(i,0:2),n_nodes)
     enddo
     if ((nods_per_elem==20) .or. (nods_per_elem==27)) then
         do i = 0,11
             ne = edges(n,i)
             ok = 1
             check_edges : do j=0,n-1
                 do k = 0,11
                     if (edges(j,k)==ne) then
                         ok = 0
                         exit check_edges
                     endif
                 enddo
             enddo check_edges
             if (ok==1) then
              if (i==0) then
                  j = 0;   k = 1
              else if (i==1) then
                  j = 1;   k = 2
              else if (i==2) then
                  j = 2;   k = 3
              else if (i==3) then
                  j = 0;   k = 3
              else if (i==4) then
                  j = 1;   k = 5
              else if (i==5) then
                  j = 4;   k = 5
              else if (i==6) then
                  j = 0;   k = 4
              else if (i==7) then
                  j = 2;   k = 6
              else if (i==8) then
                  j = 5;   k = 6
              else if (i==9) then
                  j = 6;   k = 7
              else if (i==10) then
                  j = 3;   k = 7
              else if (i==11) then
                  j = 4;   k = 7
              endif
              coord_nodes(n_nodes,:) = (coord_vertices(j,:)+coord_vertices(k,:))/2
              n_nodes = n_nodes + 1
             endif
         enddo
         if (nods_per_elem==27) then
             do i = 0,5
                 nf = faces(n,i)
                 ok = 1
                 check_faces : do j = 0,n-1
                     do k = 0,5
                         if (faces(j,k)==nf) then
                             ok = 0
                             exit check_faces
                         endif
                     enddo
                 enddo check_faces
                 if (ok==1) then
                  if (i==0) then
                      j = 0;   k = 2
                  else if (i==1) then
                      j = 0;   k = 5
                  else if (i==2) then
                      j = 1;   k = 6
                  else if (i==3) then
                      j = 2;   k = 7
                  else if (i==4) then
                      j = 0;   k = 7
                  else if (i==5) then
                      j = 4;   k = 6
                  endif
                  coord_nodes(n_nodes,:) = (coord_vertices(j,:)+coord_vertices(k,:))/2.
                  n_nodes = n_nodes + 1
                 endif
             enddo
             coord_nodes(n_nodes,:) = (coord_vertices(0,:)+coord_vertices(6,:))/2.
             n_nodes = n_nodes + 1
         endif
     endif
 enddo
 allocate (tmp2D(0:n_nodes-1,0:2))
 tmp2D(0:n_nodes-1,0:2) = coord_nodes(0:n_nodes-1,0:2)
 deallocate (coord_nodes)
 allocate (coord_nodes(0:n_nodes-1,0:2))
 coord_nodes(0:n_nodes-1,0:2) = tmp2D(0:n_nodes-1,0:2)
 deallocate (tmp2D)

 !!! Sorting the nodes !!!
 n = n_nodes-1
 ok = 1
 epsil_xi = dx/1000
 epsil_eta = dy/1000
 if (curve) then
     epsil_r = dR(n_layers-1)/1000
 else
     epsil_r = dz/1000
 endif
 do while (ok==1)
     ok = 0
     do i = 0,n-1
         if (coord_nodes(i,2)>coord_nodes(i+1,2)+epsil_r) then
             tmp1D(0:2) = coord_nodes(i+1,0:2)
             coord_nodes(i+1,0:2) = coord_nodes(i,0:2)
             coord_nodes(i,0:2) = tmp1D(0:2)
             ok = 1
         else if (egal(coord_nodes(i,2),coord_nodes(i+1,2),epsil_r)==.true.) then
             if (coord_nodes(i,1)>coord_nodes(i+1,1)+epsil_eta) then
                 tmp1D(0:2) = coord_nodes(i+1,0:2)
                 coord_nodes(i+1,0:2) = coord_nodes(i,0:2)
                 coord_nodes(i,0:2) = tmp1D(0:2)
                 ok = 1
             else if (egal(coord_nodes(i,1),coord_nodes(i+1,1),epsil_eta)==.true.) then
                     if (coord_nodes(i,0)>coord_nodes(i+1,0)+epsil_xi) then
                     tmp1D(0:2) = coord_nodes(i+1,0:2)
                     coord_nodes(i+1,0:2) = coord_nodes(i,0:2)
                     coord_nodes(i,0:2) = tmp1D(0:2)
                     ok = 1
                 else if (egal(coord_nodes(i,0),coord_nodes(i+1,0),epsil_xi)==.true.) then
                     print *,"A node appears twice !!!"
                     stop
                 endif
             endif
         endif
     enddo
     n = n-1
 enddo

 !!! Associating to each element 27 (or 20 depending on the inputs) nodes !!!
 !!! If nods_per_elem==8, the nodes are the vertices !!!
 if ((nods_per_elem==20) .or. (nods_per_elem==27)) then
     allocate (nodes(0:nelem_in_proc(proc)-1,0:nods_per_elem-1))
     do nel = 0,nelem_in_proc(proc)-1
         do i = 0,7
             nv = elmnts_local(nel,i)
             num = which_vertices(nv)
             coord_vertices(i,0) = xco(num)
             coord_vertices(i,1) = yco(num)
             coord_vertices(i,2) = zco(num)
         enddo
         dist = abs(coord_vertices(0,0) - coord_vertices(1,0))
         epsil_xi = dist / 1000
         dist = abs(coord_vertices(0,1) - coord_vertices(3,1))
         epsil_eta = dist / 1000
         dist = abs(coord_vertices(0,2) - coord_vertices(4,2))
         epsil_r = dist / 1000
         do i = 0,7
             nv = elmnts_local(nel,i)
             num = which_vertices(nv)
             search : do n = 0,n_nodes-1
                 if ( egal(coord_nodes(n,0),xco(num),epsil_xi) .and. &
                      egal(coord_nodes(n,1),yco(num),epsil_eta) .and. &
                      egal(coord_nodes(n,2),zco(num),epsil_r) ) then
                     nodes(nel,i) = n
                     exit search
                 endif
             enddo search
         enddo
         nodes(nel,8) = nodes(nel,0) + 1
         nodes(nel,10) = nodes(nel,3) + 1
         nodes(nel,16) = nodes(nel,4) + 1
         nodes(nel,18) = nodes(nel,7) + 1
         tmp1D(0) = coord_vertices(0,0)
         tmp1D(1) = (coord_vertices(0,1)+coord_vertices(3,1))/2.
         tmp1D(2) = coord_vertices(0,2)
         search11 : do n = 0,n_nodes-1
             if ( egal(coord_nodes(n,0),tmp1D(0),epsil_xi) .and. &
                  egal(coord_nodes(n,1),tmp1D(1),epsil_eta) .and. &
                  egal(coord_nodes(n,2),tmp1D(2),epsil_r) ) then
                 nodes(nel,11) = n
                 exit search11
             endif
         enddo search11
         tmp1D(0) = coord_vertices(0,0)
         tmp1D(1) = coord_vertices(0,1)
         tmp1D(2) = (coord_vertices(0,2)+coord_vertices(4,2))/2.
         search12 : do n = 0,n_nodes-1
             if ( egal(coord_nodes(n,0),tmp1D(0),epsil_xi) .and. &
                  egal(coord_nodes(n,1),tmp1D(1),epsil_eta) .and. &
                  egal(coord_nodes(n,2),tmp1D(2),epsil_r) ) then
                 nodes(nel,12) = n
                 exit search12
             endif
         enddo search12
         tmp1D(0) = coord_vertices(3,0)
         tmp1D(1) = coord_vertices(3,1)
         tmp1D(2) = (coord_vertices(3,2)+coord_vertices(7,2))/2.
         search15 : do n = 0,n_nodes-1
             if ( egal(coord_nodes(n,0),tmp1D(0),epsil_xi) .and. &
                  egal(coord_nodes(n,1),tmp1D(1),epsil_eta) .and. &
                  egal(coord_nodes(n,2),tmp1D(2),epsil_r) ) then
                 nodes(nel,15) = n
                 exit search15
             endif
         enddo search15
         tmp1D(0) = coord_vertices(4,0)
         tmp1D(1) = (coord_vertices(4,1)+coord_vertices(7,1))/2.
         tmp1D(2) = coord_vertices(4,2)
         search19 : do n = 0,n_nodes-1
             if ( egal(coord_nodes(n,0),tmp1D(0),epsil_xi) .and. &
                  egal(coord_nodes(n,1),tmp1D(1),epsil_eta) .and. &
                  egal(coord_nodes(n,2),tmp1D(2),epsil_r) ) then
                 nodes(nel,19) = n
                 exit search19
             endif
         enddo search19
         if (nods_per_elem==20) then
             nodes(nel,9) = nodes(nel,11) + 1
             nodes(nel,13) = nodes(nel,12) + 1
             nodes(nel,14) = nodes(nel,15) + 1
             nodes(nel,17) = nodes(nel,19) + 1
         else
             nodes(nel,20) = nodes(nel,11) + 1
             nodes(nel,9) = nodes(nel,11) + 2
             nodes(nel,21) = nodes(nel,12) + 1
             nodes(nel,13) = nodes(nel,12) + 2
             nodes(nel,23) = nodes(nel,15) + 1
             nodes(nel,14) = nodes(nel,15) + 2
             nodes(nel,25) = nodes(nel,19) + 1
             nodes(nel,17) = nodes(nel,19) + 2
             tmp1D(0) = coord_vertices(0,0)
             tmp1D(1) = (coord_vertices(0,1)+coord_vertices(3,1))/2.
             tmp1D(2) = (coord_vertices(0,2)+coord_vertices(4,2))/2.
             search24 : do n = 0,n_nodes-1
                 if ( egal(coord_nodes(n,0),tmp1D(0),epsil_xi) .and. &
                      egal(coord_nodes(n,1),tmp1D(1),epsil_eta) .and. &
                      egal(coord_nodes(n,2),tmp1D(2),epsil_r) ) then
                     nodes(nel,24) = n
                     exit search24
                 endif
             enddo search24
             nodes(nel,26) = nodes(nel,24) + 1
             nodes(nel,22) = nodes(nel,24) + 2
         endif
     enddo
 endif

 !!! Calculating the cartesian coordinates of each node !!!
 if (curve) then
     if (ellipticity)   allocate(Rsph(0:n_nodes-1))
     do n = 0,n_nodes-1
         R = coord_nodes(n,2)
         Y = tan (pi*coord_nodes(n,1)/180.d0)
         X = tan (pi*coord_nodes(n,0)/180.d0)
         D = sqrt(1.d0 + Y**2 + X**2)
         if (ellipticity) then
             ! Passage en coordonnees cartesiennes
             xa = X/D;   ya = Y/D;   za = 1/D
             ! Faire la rotation du chunk de ref vers le chunk reel
             xs = rot(0,0)*xa + rot(0,1)*ya + rot(0,2)*za
             ys = rot(1,0)*xa + rot(1,1)*ya + rot(1,2)*za
             zs = rot(2,0)*xa + rot(2,1)*ya + rot(2,2)*za
             ! Passage en spherique
             call cart2sph(xs,ys,zs,alpha,theta,phi)
             ! Calcul des rayons elliptiques
             P2 = (3*cos(theta)**2 - 1) / 2
             alpha = get_ellipticity(radius(n_layers-2))
             rad_loc(0) = radius(n_layers-2) * (1.d0 - 2.d0/3.d0 * alpha * P2)
             alpha = get_ellipticity(radius(n_layers-1))
             rad_loc(1) = radius(n_layers-1) * (1.d0 - 2.d0/3.d0 * alpha * P2)
             alpha = get_ellipticity(Rterre)
             rad_loc(2) = Rterre * (1.d0 - 2.d0/3.d0 * alpha * P2)
             alpha = get_ellipticity(R)
             rad_loc(3) = R * (1.d0 - 2.d0/3.d0 * alpha * P2)
         endif
         ! When the Moho and/or the surface topography are required, we change R for some nodes
         if ((model==2) .or. (model==3)) then
             if (R>limit_inf) then
                 epsil_xi = dx/1000
                 dist = dx/2
                 i = 0
                 do while (egal(coord_nodes(n,0),x_len+i*dist,epsil_xi)==.false.)
                     i = i+1
                 enddo
                 epsil_eta = dy/1000
                 dist = dy/2
                 j = 0
                 do while (egal(coord_nodes(n,1),y_len+j*dist,epsil_eta)==.false.)
                     j = j+1
                 enddo
                 if (R<limit_med) then
                     epsil_r = dR(n_layers-2)/1000
                     dist = dR(n_layers-2)/2
                     k = 1
                     do while (egal(R,radius(n_layers-2)+k*dist,epsil_r)==.false.)
                         k = k+1
                     enddo
                     dist = (Rterre-radius(n_layers-2)-depth_moho(i,j)) / (2*nz(n_layers-2))
                     R = radius(n_layers-2) + k*dist
                     if (ellipticity) then
                         Rsph(n) = R
                         dist = (rad_loc(2)-rad_loc(0)-depth_moho(i,j)) / (2*nz(n_layers-2))
                         R = rad_loc(0) + k*dist
                     endif
                     if (2*dist<dzmin(0,proc))   dzmin(0,proc) = 2*dist
                     if (2*dist>dzmax(0,proc))   dzmax(0,proc) = 2*dist
                 else
                     epsil_r = dR(n_layers-1)/1000
                     dist = dR(n_layers-1)/2
                     k = 0
                     do while (egal(R,radius(n_layers-1)+k*dist,epsil_r)==.false.)
                         k = k+1
                     enddo
                     if (topo_log) then
                         dist = (depth_moho(i,j) + altitude_topo(i,j)) / (2*nz(n_layers-1))
                     else
                         dist = depth_moho(i,j) / (2*nz(n_layers-1))
                     endif
                     R = Rterre-depth_moho(i,j) + k*dist
                     if (ellipticity) then
                         Rsph(n) = R
                         R = rad_loc(2)-depth_moho(i,j) + k*dist
                     endif
                     if (2*dist<dzmin(1,proc))   dzmin(1,proc) = 2*dist
                     if (2*dist>dzmax(1,proc))   dzmax(1,proc) = 2*dist
                 endif
             else
                 if (ellipticity) then
                     Rsph(n) = R
                     R = rad_loc(3)
                 endif
             endif
         else
             if (topo_log) then
                 if (R>limit_med) then
                     epsil_xi = dx/1000
                     dist = dx/2
                     i = 0
                     do while (egal(coord_nodes(n,0),x_len+i*dist,epsil_xi)==.false.)
                         i = i+1
                     enddo
                     epsil_eta = dy/1000
                     dist = dy/2
                     j = 0
                     do while (egal(coord_nodes(n,1),y_len+j*dist,epsil_eta)==.false.)
                         j = j+1
                     enddo
                     epsil_r = dR(n_layers-1)/1000
                     dist = dR(n_layers-1)/2
                     k = 0
                     do while (egal(R,radius(n_layers-1)+k*dist,epsil_r)==.false.)
                         k = k+1
                     enddo
                     dist = (Rterre-radius(n_layers-1)+altitude_topo(i,j)) / (2*nz(n_layers-1))
                     R = radius(n_layers-1) + k*dist
                     if (ellipticity) then
                         Rsph(n) = R
                         dist = (rad_loc(2)-rad_loc(1)+altitude_topo(i,j)) / (2*nz(n_layers-1))
                         R = rad_loc(1) + k*dist
                     endif
                 else
                     if (ellipticity) then
                         Rsph(n) = R
                         R = rad_loc(3)
                     endif
                 endif
             else
                 if (ellipticity) then
                     Rsph(n) = R
                     R = rad_loc(3)
                 endif
             endif
         endif
         coord_nodes(n,0) = R*X/D
         coord_nodes(n,1) = R*Y/D
         coord_nodes(n,2) = R/D
     enddo
 endif

 !!! Writing the meshfile !!!
 write (meshfilename(11:13), '(i3.3)') proc
 open (11, file = trim(meshfilename))
 write (11,"(1i3)") n_dim
 write (11,*)
 write (11,"(l3,i3,l3)") curve, model, ellipticity
 if (curve) then
    do i = 0,2
       write (11,*) (rot(i,j),j=0,2)
    enddo
 endif
 write (11,*)
 write (11,"(1i6)") n_nodes
 do n = 0,n_nodes-1
     write (11,*) coord_nodes(n,:)
     if (ellipticity)   write (11,*) Rsph(n)
 enddo 
 write (11,*)
 write (11,"(1i6)") nelem_in_proc(proc)
 write (11,*)
 do n = 0,nelem_in_proc(proc)-1
     nel = which_elem_in_proc(proc,n)
     write (11,"(4i6)") Material(nel), moho_position(nel), topo_position(nel), t_reversal(nel)
 enddo
 write (11,*)
 write (11,"(1i6)") nods_per_elem
 if (nods_per_elem==20) then
     do n = 0,nelem_in_proc(proc)-1
         write (11,"(20i7)") (nodes(n,i),i=0,nods_per_elem-1)
     enddo
     if (random_trac==.false. .and. output_chunk==.false.)   deallocate (nodes)
 else if (nods_per_elem==27) then
     do n = 0,nelem_in_proc(proc)-1
         write (11,"(27i7)") (nodes(n,i),i=0,nods_per_elem-1)
     enddo
     if (random_trac==.false. .and. output_chunk==.false.)  deallocate (nodes)
 else if (nods_per_elem==8) then
     do n = 0,nelem_in_proc(proc)-1
         write (11,"(8i7)") (elmnts_local(n,i),i=0,nods_per_elem-1)
     enddo
 else
     print *,"Bad number of nodes"
     stop
 endif
 write (11,*)
 write (11,"(1i6)") n_faces
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(6i6)") (faces(n,i),i=0,5)
 enddo
 write (11,*)
 write (11,"(1i6)") n_edges
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(12i6)") (edges(n,i),i=0,11)
 enddo
 write (11,*)
 write (11,"(1i6)") n_vertices
 do n = 0,nelem_in_proc(proc)-1
     write (11,"(8i6)") (elmnts_local(n,i),i=0,7)
 enddo
 write (11,*)
 write (11,"(1i6)") nparts
 do n = 0,nparts-1
     write (11,"(3i6)") nf_shared(n), ne_shared(n), nv_shared(n)
     do nf = 0,nf_shared(n)-1
         write (11,"(1i6)") faces_shared(n,nf)
     enddo
     do ne = 0,ne_shared(n)-1
         write (11,"(1i6)") edges_shared(n,ne)
     enddo
     do nv = 0,nv_shared(n)-1
         write (11,"(1i6)") vertices_shared(n,nv)
     enddo
 enddo
 close (11)

 !!! Writing chunk.out !!!
 if (output_chunk) then
     if (proc==0) then
         open(30,file="chunk_0.out")
         open(31,file="chunk_1.out")
         open(32,file="chunk_2.out")
         open(33,file="chunk_3.out")
         open(34,file="chunk_4.out")
         open(35,file="chunk_5.out")
     endif
     do k = 0,nz_tot-1
      do j = 0,ny-1
       do i = 0,nx-1
           if (k==0) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1) then
                   if (nods_per_elem==8) then
                       call write_mesh (nods_per_elem,coord_nodes(elmnts_local(n_out,:),:),0)
                   else
                       call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),0)
                   endif
               endif
           endif
           if (k==nz_tot-1) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1) then
                   if (nods_per_elem==8) then
                       call write_mesh (nods_per_elem,coord_nodes(elmnts_local(n_out,:),:),5)
                   else
                       call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),5)
                   endif
               endif
           endif
           if (j==0) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1) then
                   if (nods_per_elem==8) then
                       call write_mesh (nods_per_elem,coord_nodes(elmnts_local(n_out,:),:),1)
                   else
                       call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),1)
                   endif
               endif
           endif
           if (j==ny-1) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1) then
                   if (nods_per_elem==8) then
                       call write_mesh (nods_per_elem,coord_nodes(elmnts_local(n_out,:),:),3)
                   else
                       call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),3)
                   endif
               endif
           endif
           if (i==0) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1) then
                   if (nods_per_elem==8) then
                       call write_mesh (nods_per_elem,coord_nodes(elmnts_local(n_out,:),:),4)
                   else
                       call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),4)
                   endif
               endif
           endif
           if (i==nx-1) then
               i_count = i + j*nx + k*n_elem_xy
               n_out = -1
               call find_numloc (nelem_in_proc(proc),which_elem_in_proc(proc,:),i_count,n_out)
               if (n_out>-1) then
                   if (nods_per_elem==8) then
                       call write_mesh (nods_per_elem,coord_nodes(elmnts_local(n_out,:),:),2)
                   else
                       call write_mesh (nods_per_elem,coord_nodes(nodes(n_out,:),:),2)
                   endif
               endif
           endif
       enddo
      enddo
     enddo
     if (nods_per_elem==20 .or. nods_per_elem==27)   deallocate (nodes)
 endif


 !!! Writing the file for random forces !!!
 if (random_trac) then
     do n = 0,nelem_in_proc(proc)-1
         nel = which_elem_in_proc(proc,n)
         if (ellipticity) then
             call table_random_t (proc, n, Material(nel), coord_nodes(nodes(n,25),:), &
                                  moho_position(nel), topo_position(nel), &
                                  rot, Rterre, pi, ellipticity, Rsph(nodes(n,25)))
         else
             call table_random_t (proc, n, Material(nel), coord_nodes(nodes(n,25),:), &
                                  moho_position(nel), topo_position(nel), &
                                  rot, Rterre, pi, ellipticity)
         endif
     enddo 
     deallocate (nodes)
 endif

 !!! Deallocation !!!
 deallocate (which_vertices, elmnts_local, faces, edges, coord_nodes)
 deallocate (nf_shared, ne_shared, nv_shared, faces_shared, edges_shared, vertices_shared)
 if (ellipticity)   deallocate(Rsph)
enddo


!!! Closing the file for random forces !!!
if (random_trac)   close (25)
!!! Closing chunk_?.out !!!
if (output_chunk) then
    close (30)
    close (31)
    close (32)
    close (33)
    close (34)
    close (35)
endif


!!! Returning the smallest and the largest vertical size among the elements above and below the Moho !!!
if ((model==2) .or. (model==3)) then
    do proc = 1,nparts-1
        if (dzmin(0,proc)<dzmin(0,0))   dzmin(0,0)=dzmin(0,proc)
        if (dzmin(1,proc)<dzmin(1,0))   dzmin(1,0)=dzmin(1,proc)
        if (dzmax(0,proc)>dzmax(0,0))   dzmax(0,0)=dzmax(0,proc)
        if (dzmax(1,proc)>dzmax(1,0))   dzmax(1,0)=dzmax(1,proc)
    enddo
    write(*,"(a,f5.1,a)"),"The smallest vertical size in the crust is",dzmin(1,0)/1000," km."
    write(*,"(a,f6.1,a)"),"The smallest vertical size in the very upper mantle is",dzmin(0,0)/1000," km."
    if (dzmax(1,0)>dRmax) then
        write(*,"(a,f5.1,a)"),"WARNING: The largest element in the crust is larger than &
                               the horizontal step (dzmax=",dzmax(1,0)/1000," km)."
    endif
    if (dzmax(0,0)>dRmax*vs(nb_couches_model-2)/vs(nb_couches_model-1)) then
        write(*,"(a,f6.1,a)"),"WARNING: The largest element in the very upper mantle is larger than &
                               the horizontal step (dzmax=",dzmax(0,0)/1000," km)."
    endif
    deallocate (dzmin, dzmax, depth_moho)
endif


!!! Deallocation !!!
do nel = 0,n_elem-1
    if (part(nel) /= nparts-1) then
        do proc = part(nel)+1,nparts-1
            deallocate (memory(nel)%rank(proc)%Obj)
        enddo
        deallocate (memory(nel)%rank)
    endif
enddo
deallocate (xco, yco, zco, moho_position, topo_position)
deallocate (elmnts, adjwgt, dxadj, dxadjncy, vwgt, part)
deallocate (Material, nelem_in_proc, which_elem_in_proc, memory)
if (curve) deallocate (radius, nz, dR)


contains


 logical function egal(a,b,epsil)

 doubleprecision :: a,b,epsil

 egal = .false.
 if (abs(a-b) < epsil) egal = .true.

 end function


 doubleprecision function deg2rad(val)

 doubleprecision :: val
 doubleprecision, parameter :: pi = 3.141592653

 deg2rad = pi*val/180.

 end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bottom_remarks(rayon,rad,i,dx,Rterre,model)

implicit none

doubleprecision, intent(IN) :: dx, Rterre
doubleprecision, intent(INOUT) :: rayon
doubleprecision, dimension(0:11), intent(IN) :: rad
integer, intent(IN) :: i, model
doubleprecision :: dr
doubleprecision, parameter :: pi = 3.14159265
integer :: ok


dr = pi*dx*Rterre/180
ok = 0
if (i<2) then
    print *,"You are under the CMB. The elements of the core will be PML."
    rayon = rad(1)-dr
    ok = 1
else if (rad(i)-rayon<dr) then
    if (rad(i)-dr<rad(i-1)) then
        rayon = rad(i-1)
    else
        rayon = rad(i)-dr
    endif
    ok = 1
endif

if ((model==2 .or. model==3) .and. (rayon>6291000)) then
    print *,"If you want the Moho, the bottom has to be deeper than 80 km."
    stop
endif

if (ok==1) write(*,"(a,f8.0)"),"Optimization of the bottom: depth = ", Rterre - rayon

end subroutine bottom_remarks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sort(vector, length)

implicit none

integer, intent(IN) :: length
integer, dimension(0:length-1), intent(INOUT) :: vector
integer :: n, ok, i, tmp


n = length-1
ok = 1
do while (ok==1)
    ok = 0
    do i = 0,n-1
        if (vector(i)>vector(i+1)) then
            tmp = vector(i+1)
            vector(i+1) = vector(i)
            vector(i) = tmp
            ok = 1
        endif
    enddo
    n = n-1
enddo

return
end subroutine sort


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_vertices(nv,n,elmnts_local,coord,nod,nb)

implicit none

integer, intent(IN) :: n ,nv
integer, intent(INOUT) :: nb
integer, dimension(0:n-1,0:7), intent(IN) :: elmnts_local
doubleprecision, dimension(0:2), intent(IN) :: nod
doubleprecision, dimension(0:2), intent(INOUT) :: coord

integer :: i, j, ok


ok = 1
loop : do i = 0,n-1
    do j = 0,7
        if (elmnts_local(i,j)==nv) then
            ok = 0
            exit loop
        endif
    enddo
enddo loop
if (ok==1) then
    coord(0:2) = nod(0:2)
    nb = nb+1
endif

end subroutine check_vertices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine define_model(num,nb_couches_model,rad,vs)

implicit none

integer, intent(IN) :: num
integer, intent(OUT) :: nb_couches_model
doubleprecision, dimension(0:11), intent(OUT) :: rad, vs


rad = -1;   vs = -1

if (num==1) then

    nb_couches_model = 12

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
    rad(8)=6291000; rad(9)=6346600; rad(10)=6356000

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4643.90
    vs(8)=4470.52; vs(9)=4491.01; vs(10)=3900.00

else if (num==2) then

    nb_couches_model = 10

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
    rad(8)=6346600

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4643.90
    vs(8)=4491.01

else if (num==3) then

    nb_couches_model = 9

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6346600

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4491.01

else if (num==5) then

    nb_couches_model = 7

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5701000
    rad(4)=5971000; rad(5)=6311000

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=5920.61
    vs(4)=4877.09; vs(5)=4556.55

else if (num==6) then

    nb_couches_model = 11

    rad(0)=1221500; rad(1)=3480000; rad(2)=3630000; rad(3)=5600000
    rad(4)=5701000; rad(5)=5771000; rad(6)=5971000; rad(7)=6151000
    rad(8)=6291000; rad(9)=6346600

    vs(0)=3054.31; vs(1)=0.00; vs(2)=7265.97; vs(3)=6240.39
    vs(4)=5945.13; vs(5)=5516.02; vs(6)=4932.49; vs(7)=4643.90
    vs(8)=4470.52; vs(9)=4491.01

else

    print *,"ERROR: ",num, " is a wrong input"
    stop

endif

end subroutine define_model


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_moho(theta_rad,phi_rad,moho,m,n,latmin,latmax,dlat,lonmin,lonmax,dlong,depth)

implicit none

integer, intent(IN) :: m,n 
doubleprecision, intent(IN) :: latmin,latmax,dlat,lonmin,lonmax,dlong
doubleprecision, intent(OUT) :: depth
doubleprecision, intent(INOUT) :: theta_rad,phi_rad
doubleprecision, dimension(0:m-1,0:n-1), intent(IN) :: moho

integer :: i,j, k,l,kk,ll, u,v,uu,vv, ok, u1,u2,v1,v2
doubleprecision :: colatmax, theta_deg,phi_deg, phi_loc, sigma, norm, epsil, Wlong,Wlat
doubleprecision, dimension(0:20) :: theta_loc, pt_inter
doubleprecision, dimension(0:20,0:20) :: pt
doubleprecision, parameter :: pi = 3.141592653


epsil = dlong/1000.d0
theta_deg = 180.d0*theta_rad/pi
phi_deg = 180.d0*phi_rad/pi

j = 0
do while (lonmin+j*dlong < phi_deg-dlong/2.d0)
    j = j + 1
enddo
if (j<2 .or. j>n-3) then
    if (dabs((lonmax+dlong-360.d0)-lonmin) > epsil) then
        write (*,'(a,f6.1)'),"PAS ASSEZ DE DONNEES POUR DETERMINER LA TOPO A LA LONGITUDE",phi_deg
        stop
    endif
endif
if (j==n) j = 0

colatmax = 90.d0 - latmin
i = 0
do while (colatmax-i*dlat > theta_deg+dlat/2.d0)
    i = i + 1
enddo
ok = 0
if (i<2 .or. i>m-3) then
    if (phi_deg >= lonmin+(n-1)*dlong) then
        v1 = n-1;   v2 = 0
        Wlong = (phi_deg - (lonmin+(n-1)*dlong)) / dlong
    else if (phi_deg<=lonmin) then
        v1 = n-1;   v2 = 0
        Wlong = (dlong - (lonmin-phi_deg)) / dlong
    else
        if (lonmin+j*dlong >= phi_deg) then
            v1 = j-1;   v2 = j
        else
            v1 = j;   v2 = j+1
        endif
        Wlong = (phi_deg - (lonmin+v1*dlong)) / dlong
    endif
    if (i==1) then
        if (colatmax-i*dlat >= theta_deg) then
            u1 = 1;   u2 = 2
        else
            u1 = 0;   u2 = 1
        endif
        Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
    else if (i==0) then
        if (colatmax-i*dlat >= theta_deg) then
            u1 = 0;   u2 = 1
            Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
        else
            if (colatmax+dlat > 180.d0) then
                u1 = 0;   u2 = 0
                Wlat = 0.d0
            else
                write (*,'(a,f5.1)'),"PAS ASSEZ DE DONNEES POUR DETERMINER LA TOPO A LA COLATITUDE",theta_deg
                stop
            endif
        endif
    else if (i==m-2) then
        if (colatmax-i*dlat <= theta_deg) then
            u1 = m-3;   u2 = m-2
        else
            u1 = m-2;   u2 = m-1
        endif
        Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
    else if (i==m-1) then
        if (colatmax-i*dlat <= theta_deg) then
            u1 = m-2;   u2 = m-1
            Wlat = ((colatmax-u1*dlat) - theta_deg) / dlat
        else
            if (colatmax-m*dlat < 0.d0) then
                u1 = m-1;   u2 = m-1
                Wlat = 0.d0
            else
                write (*,'(a,f5.1)'),"PAS ASSEZ DE DONNEES POUR DETERMINER LA TOPO A LA COLATITUDE",theta_deg
                stop
            endif
        endif
    endif
    depth = (1-Wlat) * ((1-Wlong)*moho(u1,v1) + Wlong*moho(u1,v2)) &
            + Wlat * ((1-Wlong)*moho(u2,v1) + Wlong*moho(u2,v2))
    ok = 1
endif


if (ok==0) then

 if (j-2<0) then ! j=0 ou j=1
    pt(0,0) = moho(i-2,n+(j-2))
    pt(5,0) = moho(i-1,n+(j-2))
    pt(10,0) = moho(i,n+(j-2))
    pt(15,0) = moho(i+1,n+(j-2))
    pt(20,0) = moho(i+2,n+(j-2))
 else ! j>1
    pt(0,0) = moho(i-2,j-2)
    pt(5,0) = moho(i-1,j-2)
    pt(10,0) = moho(i,j-2)
    pt(15,0) = moho(i+1,j-2)
    pt(20,0) = moho(i+2,j-2)
 endif

 if (j-1<0) then ! j=0
    pt(0,5) = moho(i-2,n-1)
    pt(5,5) = moho(i-1,n-1)
    pt(10,5) = moho(i,n-1)
    pt(15,5) = moho(i+1,n-1)
    pt(20,5) = moho(i+2,n-1)
 else ! j>0
    pt(0,5) = moho(i-2,j-1)
    pt(5,5) = moho(i-1,j-1)
    pt(10,5) = moho(i,j-1)
    pt(15,5) = moho(i+1,j-1)
    pt(20,5) = moho(i+2,j-1)
 endif

 pt(0,10) = moho(i-2,j)
 pt(5,10) = moho(i-1,j)
 pt(10,10) = moho(i,j)
 pt(15,10) = moho(i+1,j)
 pt(20,10) = moho(i+2,j)

 if (j+1>n-1) then ! j=n-1
    pt(0,15) = moho(i-2,0)
    pt(5,15) = moho(i-1,0)
    pt(10,15) = moho(i,0)
    pt(15,15) = moho(i+1,0)
    pt(20,15) = moho(i+2,0)
 else ! j<n-1
    pt(0,15) = moho(i-2,j+1)
    pt(5,15) = moho(i-1,j+1)
    pt(10,15) = moho(i,j+1)
    pt(15,15) = moho(i+1,j+1)
    pt(20,15) = moho(i+2,j+1)
 endif

 if (j+2>n-1) then ! j=n-1 ou j=n-2
    pt(0,20) = moho(i-2,(j+2)-n)
    pt(5,20) = moho(i-1,(j+2)-n)
    pt(10,20) = moho(i,(j+2)-n)
    pt(15,20) = moho(i+1,(j+2)-n)
    pt(20,20) = moho(i+2,(j+2)-n)
 else ! j<n-2
    pt(0,20) = moho(i-2,j+2)
    pt(5,20) = moho(i-1,j+2)
    pt(10,20) = moho(i,j+2)
    pt(15,20) = moho(i+1,j+2)
    pt(20,20) = moho(i+2,j+2)
 endif

 do l = 0,4
     u = 5*l
     do k = 0,4
         v = 5*k
         do ll = 0,4
             uu = u + (ll-2)
             if (uu>=0) then
                 do kk = 0,4
                     vv = v + (kk-2)
                     if (vv>=0) then
                         pt(uu,vv) = pt(u,v)
                     endif
                  enddo
             endif
         enddo
     enddo
 enddo

 sigma = dlong/2.d0
 if (j==0 .and. phi_deg>(n-1)*dlong)   j = n
 do l = 0,20
     theta_loc(l) = colatmax - (i*dlat + (l-10)*dlat/5.d0)
     do k = 0,20
         phi_loc = j*dlong + (k-10)*dlong/5.d0
         pt(l,k) = pt(l,k) * exp (- (phi_loc-phi_deg)**2 / (2.d0*sigma**2))
     enddo
 enddo
 pt_inter(:) = 0
 do l = 0,20
     do k = 0,9
         pt_inter(l) = pt_inter(l) + (pt(l,2*k) + 4*pt(l,2*k+1) + pt(l,2*k+2))
     enddo
     norm = dlong / (15*sigma*sqrt(2*pi))
     pt_inter(l) = pt_inter(l) * norm
 enddo

 sigma = dlat/2.d0
 do l = 0,20
     pt_inter(l) = pt_inter(l) * exp (- (theta_loc(l)-theta_deg)**2 / (2.d0*sigma**2))
 enddo
 depth = 0
 do k = 0,9
     depth = depth + (pt_inter(2*k) + 4*pt_inter(2*k+1) + pt_inter(2*k+2))
 enddo
 norm = dlat / (15*sigma*sqrt(2*pi))
 depth = depth * norm

endif

end subroutine read_moho


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program mesher
