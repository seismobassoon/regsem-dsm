subroutine one_more_predic (Tdomain, rg, ntime)


use sdomains

implicit none

type (domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: rg, ntime

logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
integer :: n, mat, elem, nf,ne,nv, i, i_simu
doubleprecision :: dt, alpha,bega,gam1


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

call rec_mirror (Tdomain, ntime, rg)


return
end subroutine one_more_predic
