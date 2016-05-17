subroutine rec_mirror (Tdomain, ntime, rg)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, rg

integer :: len_DP, n, i, x,y, comp
doubleprecision, dimension(0:Tdomain%recl_mirror-1) :: tmp
character*60 :: mirror_file


if (ntime==1) then
    inquire(iolength=len_DP) tmp
    if (len_DP/=0) then
     write (mirror_file,"(a,I3.3)") "mirror",rg
     open (44,file=mirror_file,access='direct',form="unformatted",recl=len_DP,status='replace')
    endif
endif

i = 0
do n = 0, Tdomain%n_face-1
    if (Tdomain%sFace(n)%mirror) then
        do x = 1, Tdomain%sFace(n)%ngll1-2
         do y = 1, Tdomain%sFace(n)%ngll2-2
          do comp = 0,2
              tmp(i) = Tdomain%sFace(n)%sSimu(0)%Forces(x,y,comp)
              i = i + 1
          enddo
         enddo
        enddo
    endif
enddo
do n = 0,Tdomain%n_edge-1
    if (Tdomain%sEdge(n)%mirror) then
        do x = 1, Tdomain%sEdge(n)%ngll-2
         do comp = 0,2
             tmp(i) = Tdomain%sEdge(n)%sSimu(0)%Forces(x,comp)
             i = i + 1
         enddo
        enddo
    endif
enddo
do n = 0,Tdomain%n_vertex-1
    if (Tdomain%sVertex(n)%mirror) then
        do comp = 0,2
            tmp(i) = Tdomain%sVertex(n)%sSimu(0)%Forces(comp)
            i = i + 1
        enddo
    endif
enddo
if (i/=Tdomain%recl_mirror)   stop 'A TMP BUFFER HAS A LENGTH DIFFERENT FROM RECL_MIRROR !!!'
if (Tdomain%recl_mirror/=0)   write (44,rec=ntime) tmp

if (ntime==Tdomain%sTimeParam%ntime .and. Tdomain%recl_mirror/=0)   close (44)


end subroutine rec_mirror
