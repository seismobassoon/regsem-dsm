subroutine cart2sph (x,y,z,r,theta,phi)

implicit none

doubleprecision, intent(IN) :: x,y,z
doubleprecision, intent(OUT) :: r,theta,phi
doubleprecision :: dx
doubleprecision, parameter :: pi = 3.141592653

r = dsqrt(x**2 + y**2 + z**2)
if (r==0) then
    theta = 0;   phi = 0
else
    dx = z/r
    if (dx >= 1) then
        theta = 0
    else if (dx <= -1) then
        theta = pi
    else
        theta = acos(dx)
    endif
    if ((theta==0) .or. (theta==pi)) then
        phi = 0
    else
        dx = x/(r*dsin(theta))
        if (dx > 1) then
            phi = 0
        else if (dx < -1) then
            phi = pi
        else
            phi = acos(dx)
            if (y < 0)   phi = 2*pi - phi
        endif
    endif
endif

return
end subroutine cart2sph
