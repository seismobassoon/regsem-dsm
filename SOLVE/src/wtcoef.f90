!-----------------------------------------------------
 subroutine wtcoef(f,f1,f2,f3,f4,wt)
!-----------------------------------------------------
   implicit none

   doubleprecision, intent(in) ::  f,f1,f2,f3,f4
   doubleprecision, intent(out)::  wt
   doubleprecision, parameter :: pi = 3.141592653

   if (f3.gt.f4) stop 'wtcoef: f3>f4 '
   if (f1.gt.f2) stop 'wtcoef: f1>f2 '
   if (f.le.f3.and.f.ge.f2) then
      wt = 1.
   else if (f.gt.f4.or.f.lt.f1 ) then
      wt = 0.
   else if (f.gt.f3.and.f.le.f4) then
      wt = 0.5*(1.0+cos(pi*(f-f3)/(f4-f3)))
   else if (f.ge.f1.and.f.lt.f2) then
      wt = 0.5*(1.0+cos(pi*(f-f2)/(f2-f1)))
   endif
!-----------------------------------------------------
 end subroutine wtcoef
!-----------------------------------------------------
