subroutine mirror_definition (Tdomain)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain

integer :: n, i, nf,ne,nv, x,y


do n = 0,Tdomain%n_face-1
   Tdomain%sFace(n)%mirror = .false.
enddo
do n = 0,Tdomain%n_edge-1
   Tdomain%sEdge(n)%mirror = .false.
enddo
do n = 0,Tdomain%n_vertex-1
   Tdomain%sVertex(n)%mirror = .false.
enddo

do n = 0,Tdomain%n_elem-1
 i = Tdomain%specel(n)%mirror_position
 if (i==1 .or. i==7 .or. i==8 .or. i==9 .or. i==10 .or. i==19 .or. i==20 .or. i==21 .or. i==22) then
    nf = Tdomain%specel(n)%Near_Faces(0);   Tdomain%sFace(nf)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%mirror = .true.
 endif
 if (i==2 .or. i==7 .or. i==11 .or. i==12 .or. i==15 .or. i==19 .or. i==20 .or. i==23 .or. i==24) then
    nf = Tdomain%specel(n)%Near_Faces(1);   Tdomain%sFace(nf)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%mirror = .true.
 endif
 if (i==3 .or. i==8 .or. i==12 .or. i==13 .or. i==16 .or. i==20 .or. i==21 .or. i==24 .or. i==25) then
    nf = Tdomain%specel(n)%Near_Faces(2);   Tdomain%sFace(nf)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%mirror = .true.
 endif
 if (i==4 .or. i==9 .or. i==13 .or. i==14 .or. i==17 .or. i==21 .or. i==22 .or. i==25 .or. i==26) then
    nf = Tdomain%specel(n)%Near_Faces(3);   Tdomain%sFace(nf)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%mirror = .true.
 endif
 if (i==5 .or. i==10 .or. i==11 .or. i==14 .or. i==18 .or. i==19 .or. i==22 .or. i==23 .or. i==26) then
    nf = Tdomain%specel(n)%Near_Faces(4);   Tdomain%sFace(nf)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%mirror = .true.
 endif
 if (i==6 .or. i==15 .or. i==16 .or. i==17 .or. i==18 .or. i==23 .or. i==24 .or. i==25 .or. i==26) then
    nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%mirror = .true.
    ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%mirror = .true.
    nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%mirror = .true.
 endif
enddo

i = 0
do n = 0,Tdomain%n_face-1
   if (Tdomain%sFace(n)%mirror) then
      do x = 1, Tdomain%sFace(n)%ngll1-2
         do y = 1, Tdomain%sFace(n)%ngll2-2
            i = i + 1
         enddo
      enddo
   endif
enddo
do n = 0,Tdomain%n_edge-1
   if (Tdomain%sEdge(n)%mirror) then
      do x = 1, Tdomain%sEdge(n)%ngll-2
         i = i + 1 
      enddo
   endif
enddo
do n = 0,Tdomain%n_vertex-1
   if (Tdomain%sVertex(n)%mirror)   i = i + 1
enddo
Tdomain%recl_mirror = 3*i


return
end subroutine mirror_definition

