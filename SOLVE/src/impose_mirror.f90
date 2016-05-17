subroutine impose_mirror (Tdomain, ntime, rg)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: ntime, rg

integer :: len_DP, n, i, x,y, comp
doubleprecision :: current_t, t_max, wtaper
doubleprecision, dimension(0:Tdomain%recl_mirror-1) :: tmp
character*60 :: mirror_file

doubleprecision, parameter :: taperwidth = 200.d0   ! Mettre 0.d0 pour ne pas avoir de taper


if (ntime==0) then
    inquire(iolength=len_DP) tmp
    if (len_DP/=0) then
     write (mirror_file,"(a,I3.3)") "mirror",rg
     open (44,file=mirror_file,access='direct',form="unformatted",recl=len_DP,status='old')
    endif
endif

if (Tdomain%recl_mirror/=0)   read(44,rec=Tdomain%sTimeParam%ntime-ntime) tmp

! Taper applique au debut du mirroir pour eviter une arrivee trop brusque de l'energie.
! Ca peut arriver si le mirroir n'est pas bien a zero.
current_t = Tdomain%sTimeParam%rtime
t_max = Tdomain%sTimeParam%duration + Tdomain%sTimeParam%dt   ! On ajoute un dt pour etre bien sur que t est trjs < t_max
call wtcoef(current_t,0.d0,taperwidth,t_max,t_max,wtaper)
tmp = wtaper*tmp

i = 0
do n = 0, Tdomain%n_face-1
    if (Tdomain%sFace(n)%mirror) then
        do x = 1, Tdomain%sFace(n)%ngll1-2
         do y = 1, Tdomain%sFace(n)%ngll2-2
          do comp = 0,2
              Tdomain%sFace(n)%sSimu(0)%Forces(x,y,comp) = tmp(i)
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
             Tdomain%sEdge(n)%sSimu(0)%Forces(x,comp) = tmp(i)
             i = i + 1
         enddo
        enddo
    endif
enddo
do n = 0,Tdomain%n_vertex-1
    if (Tdomain%sVertex(n)%mirror) then
        do comp = 0,2
            Tdomain%sVertex(n)%sSimu(0)%Forces(comp) = tmp(i)
            i = i + 1
        enddo
    endif
enddo
if (i/=Tdomain%recl_mirror)   stop 'A TMP BUFFER HAS A LENGTH DIFFERENT FROM RECL_MIRROR !!!'

if (ntime==Tdomain%sTimeParam%ntime-1 .and. Tdomain%recl_mirror/=0)   close (44)


end subroutine impose_mirror
