module DerivatorMOD
  use tools
  implicit none
  private
  public :: getDeriv
  integer :: VSIZE = 3
  procedure(derivBackward), pointer :: pMethod

  interface getDeriv
     module procedure getDeriv
  end interface getDeriv
contains
  real(dp) function getDeriv(vect, h, method)
    implicit none
    real(dp), dimension(VSIZE), intent(in) :: vect
    real(dp) :: h
    character(*) :: method
    call chooseMethod(method)
    getDeriv = pMethod(vect, h)
  end function getDeriv

  subroutine chooseMethod(method)
    implicit none
    character(*) :: method
    select case(method)
    case('b')
       pMethod => derivBackward
    case('backward')
       pMethod => derivBackward
    case('c')
       pMethod => derivCentered
    case('centered')
       pMethod => derivCentered
    case('f')
       pMethod => derivForward
    case('forward')
       pMethod => derivForward
    case default
    end select
  end subroutine chooseMethod

  real(dp) function derivBackward(vect, h)
    implicit none
    real(dp), dimension(VSIZE), intent(in) :: vect
    real(dp) :: h
    derivBackward = (3*vect(3)-4*vect(2)+vect(1))/(2*h)
  end function derivBackward
  real(dp) function derivCentered(vect, h)
    implicit none
    real(dp), dimension(VSIZE), intent(in) :: vect
    real(dp) :: h
    derivCentered = (vect(3)-vect(1))/(2*h)
  end function derivCentered
  real(dp) function derivForward(vect, h)
    implicit none
    real(dp), dimension(VSIZE), intent(in) :: vect
    real(dp) :: h
    derivForward = (-vect(3)+4*vect(2)-3*vect(1))/(2*h)
  end function derivForward
end module DerivatorMOD
  
