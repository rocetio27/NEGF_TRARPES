function r3cross(x,y)
implicit none
real(8) :: x(3),y(3)
real(8) :: r3cross(3)
r3cross(1)=x(2)*y(3)-x(3)*y(2)
r3cross(2)=x(3)*y(1)-x(1)*y(3)
r3cross(3)=x(1)*y(2)-x(2)*y(1)
end function
