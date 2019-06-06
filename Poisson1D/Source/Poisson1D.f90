Module Poisson1DMod
  Use tools
  Use domainMod
  Use element1DMod
  Use GaussLegendreFEM
  Use sparseMod
  Use Print_All
  Implicit none
  Private
  Public :: Poisson1DType

  Type Poisson1DType
     Type(DomainType) :: domain
     Procedure(fDefault), pointer, nopass :: kCoef
     Procedure(fDefault), pointer, nopass :: source
     Procedure(fDefault), pointer, nopass :: rho
     Procedure(fDefault), pointer, nopass :: c
     Type(sparseType) :: K
     Type(sparseType) :: M
     Type(sparseType) :: MLumped
     Type(sparseType) :: MInverse     
     Real(dp), Dimension(:), Allocatable :: phi
     Real(dp), Dimension(:), Allocatable :: rhs
     Real(dp), Dimension(:), Allocatable :: flux
     Integer, Dimension(:), Allocatable :: DirichletNode, NeumannNode
     Real(dp), Dimension(:), Allocatable :: DirichletValue, NeumannValue
     Real(dp) :: t0, tf, theta
     Integer :: nGauss
     Character(250) :: method, SolveMethod
   Contains
     Procedure, Public :: init
     Procedure, Public :: update
     Procedure, Public :: solve
     Procedure, Public :: solveSteadyState
     Procedure, Public :: solveTransientState
     Procedure, Private :: computeK
     Procedure, Private :: computeM
     Procedure, Private :: computeMLumped
     Procedure, Private :: getInverse
     Procedure, Private :: computeRHS
     Procedure, Private :: boundaryConditions
     Procedure, Public :: computeFlux
     Procedure, Private :: valueJacobianGaussPoints
     Procedure, Private :: valueJacobianNodePoints
     Procedure, Public :: TransientMethod
!!$     Procedure, Public :: Runge_Kutta4
  End type Poisson1DType
  
  Real(dp), Dimension(:,:,:), Allocatable :: shapeFuncMat, derivShapeFuncMat
  Real(dp), Dimension(:), Allocatable :: QVect, kVect, rhoVect, cVect
  Real(dp) :: Kij, fi, Mij
  Integer :: i, j, iElem, ip, jp
  
Contains

  Real(dp) Function fDefault(i,x)
    Implicit none
    Integer, Intent(In) :: i
    Real(dp), Intent(In) :: x
  End Function fDefault

  Subroutine init(this)
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    Call setGauss(this%nGauss)

    Allocate(shapeFuncMat(3,4,this%nGauss))
    Allocate(derivShapeFuncMat(3,4,this%nGauss))

    Call valueGauss(lineardN, derivShapeFuncMat(1,1:2,1:this%nGauss))
    Call valueGauss(cuadraticdN, derivShapeFuncMat(2,1:3,1:this%nGauss))
    Call valueGauss(cubicdN, derivShapeFuncMat(3,1:4,1:this%nGauss))
    Call valueGauss(linearN, shapeFuncMat(1,1:2,1:this%nGauss))
    Call valueGauss(cuadraticN, shapeFuncMat(2,1:3,1:this%nGauss))
    Call valueGauss(cubicN, shapeFuncMat(3,1:4,1:this%nGauss))
    Do iElem = 1, this%domain%nElem
       Call this%valueJacobianGaussPoints
    End Do
  End Subroutine init

  Subroutine update(this)
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    Call setGauss(this%nGauss)
    If(Allocated(shapeFuncMat)) then
       Deallocate(shapeFuncMat)
       Deallocate(derivShapeFuncMat)
    End If
    
    Allocate(shapeFuncMat(3,4,this%nGauss))
    Allocate(derivShapeFuncMat(3,4,this%nGauss))

    Call valueGauss(lineardN, derivShapeFuncMat(1,1:2,1:this%nGauss))
    Call valueGauss(cuadraticdN, derivShapeFuncMat(2,1:3,1:this%nGauss))
    Call valueGauss(cubicdN, derivShapeFuncMat(3,1:4,1:this%nGauss))
    Call valueGauss(linearN, shapeFuncMat(1,1:2,1:this%nGauss))
    Call valueGauss(cuadraticN, shapeFuncMat(2,1:3,1:this%nGauss))
    Call valueGauss(cubicN, shapeFuncMat(3,1:4,1:this%nGauss))
    Do iElem = 1, this%domain%nElem
       Deallocate(this%domain%element(iElem)%jacobian)
       Call this%valueJacobianGaussPoints
    End Do
  End Subroutine update

  Subroutine Solve(this)
    Implicit none
    Class(Poisson1DType) :: this
    Select Case(this%method)
    Case('Steady_State')
       Call this%solveSteadyState
    Case('Transient_State')
       Call this%solveTransientState
    End Select
  End Subroutine Solve

  Subroutine solveSteadyState(this)
    Use GMRES
    Implicit none
    Class(Poisson1DType) :: this
    Integer, parameter :: ITR_MAX = 1000
    Real(dp), parameter :: TOL_ABS = 1d-15
    Allocate(this%flux(size(this%domain%x)))
    Call this%computeK
    Call this%computeRHS
    Allocate(this%phi(this%domain%nNodes))
    Call this%boundaryConditions
    this%phi = 0.d0
!!$    Do i = 1, this%domain%nNodes
!!$       Do j = this%K%AI(i), this%K%AI(i+1)-1
!!$          Print'(A,I0,A,I0,A,E10.4)', 'k(', i, ',', this%K%AJ(j), ') = ', this%K%A(j)
!!$       End Do
!!$    End Do
!!$    Do i = 1, this%domain%nNodes
!!$       Print'(A,I0,A,E10.4)', 'f(',i,') = ', this%rhs(i)
!!$    End Do
    Call pmgmres_ilu_cr(               &
           n = this%domain%nNodes      &
         , nz_num = size(this%K%A)     &
         , ia = this%K%AI              &
         , ja = this%K%AJ              &
         , a = this%K%A                &
         , x = this%phi                &
         , rhs = this%rhs              &
         , itr_max = ITR_MAX           &
         , mr = 500                    &
         , tol_abs = TOL_ABS           &
         , tol_rel = 1d-15             &
         , verbose = .true.           )
    Call this%computeFlux
  End Subroutine solveSteadyState

  Subroutine solveTransientState(this)
    Implicit none
    Class(Poisson1DType) :: this
    Select case(this%solveMethod)
    Case('Explicit')
       this%theta = 0.
    Case('Crank-Nicolson')
       this%theta = .5
    Case('Galerkin')
       this%theta = 2./3
    Case('Implicit')
       this%theta = 1.
!!$    Case('Runge-Kutta')
!!$       Call this%Runge_Kutta4
!!$       stop
    End Select
       Call this%TransientMethod
  End Subroutine solveTransientState

  Subroutine TransientMethod(this)
    Use GMRES    
    Implicit none
    Class(Poisson1DType) :: this
    Type(sparseType) :: worksparse1
    Type(sparseType) :: worksparse2
    Character (len = 250) :: projectName, path
    Integer, Parameter :: projectData= 1, results= 2, Data= 3
    integer :: step1, step2
    Real(sp) :: time
    Real(dp) :: Delta_t, w(this%domain%nNodes)
!!$-----------------------------------------------------------------
    Open(projectData, file = 'projectData.dat')
    Read(projectData,'(250A)') projectName
    Read(projectData,'(250A)') path
    Close(projectData)
!!$-----------------------------------------------------------------
!!$--------------- Delta t---------------------------------------------
    Delta_t = (.5*this%c(this%domain%element(1)%materialIndex,1.d0)&
            *this%rho(this%domain%element(1)%materialIndex,1.d0)&
            *(abs(this%domain%x(1)-this%domain%x(2)))**2&
            /this%kcoef(this%domain%element(1)%materialIndex,1.d0))!*.3 para matriz completa
    Print*, 'Delta_t = ',Delta_t
!!$------------------------------------------------------------------
    Open(Data, file='result')
    Open(results, file = trim(projectName)//'.flavia.res')
    write(results,*)      'GiD Post Result File 1.0'
!!$------------------------------------------------------------------
    Call this%computeK
    Call worksparse1%init(nnz = size(this%K%A), rows = this%domain%nNodes)
    Call worksparse2%init(nnz = size(this%K%A), rows = this%domain%nNodes)
        Allocate(this%phi(this%domain%nNodes))
        Allocate(this%flux(size(this%domain%x)))
    Call this%computeMlumped
!se puede cambiar por m completa pero delta t tiene q ser menor
    Call aplb(                               &        
         nrow = this%domain%nNodes           &
         ,ncol = this%domain%nNodes          &
         ,a = this%Mlumped%A/Delta_t         &
         ,ja = this%Mlumped%AJ               &
         ,ia = this%Mlumped%AI               &
         ,b = this%K%A*(this%theta-1)        &
         ,jb = this%K%AJ                     &
         ,ib =this%K%AI                      &
         ,c =  worksparse1%A                 &
         ,jc = worksparse1%AJ                &
         ,ic = worksparse1%AI               )
    Call aplb(                               &       
         nrow = this%domain%nNodes           &
         ,ncol = this%domain%nNodes          &
         ,a = this%Mlumped%A/Delta_t         &
         ,ja = this%Mlumped%AJ               &
         ,ia = this%Mlumped%AI               &
         ,b = this%K%A*this%theta            &
         ,jb = this%K%AJ                     &
         ,ib = this%K%AI                     &
         ,c =  worksparse2%A                 &
         ,jc = worksparse2%AJ                &
         ,ic = worksparse2%AI               )
    Call this%computeRHS
    Call this%getInverse (this%domain%nNodes,worksparse2,this%MInverse)
    ! para m completa sacar la lumped
!!$-------------------------------------------------------------------
!!$ Initial conditions
!!$-------------------------------------------------------------------  
    this%phi = 0.d0
    time = this%t0
    step1 = 0
    step2 = 0
    Do i = 1, size(this%NeumannNode)
       this%rhs(this%NeumannNode(i)) = &
            this%rhs(this%NeumannNode(i))&
            - this%NeumannValue(i)
    End Do
    this%rhs = this%rhs(this%domain%nNodes:1:-1)
    Do i = 1, size(this%DirichletNode)
       this%phi(this%DirichletNode(i)) = this%DirichletValue(i)
    End Do
    call this%computeFlux
!!$--------------------------------------------------------------------
!!$ calculate
!!$--------------------------------------------------------------------
    Do while (time <= this%tf)
!!$------------------------------------------------------------
!!$--------------------------writing---------------------------
!!$------------------------------------------------------------
if (step1 == step2) then
       Do i = 1, this%domain%nNodes
          Write(Data,*) this%domain%x(i),this%phi(i)
       End Do
       Write(Data,*)
       Write(Data,*)
!!$------------------------------------------------------------
!!$------------------------------------------------------------
!!$------------------------------------------------------------
       write(results,*)'Result "Temperature"',' "',trim(projectName),'" ',time,'Scalar OnNodes'
       write(results,*)      'Values' 
       Do i = 1, this%domain%nNodes
          Write(results,*)      i, this%phi(i)
       End Do
       write(results,*)      'End Values'
       write(results,*)
       write(results,*)'Result "Flux"',' "',trim(projectName),'" ',time,'Scalar OnNodes'
       write(results,*)      'Values' 
       Do i = 1, this%domain%nNodes
          Write(results,*)      i, this%Flux(i)
       End Do
       write(results,*)      'End Values'
       write(results,*)
       step2 = step2+10
end if
!!$------------------------------------------------------------
!!$----------------------end writing---------------------------
!!$------------------------------------------------------------
       time = time + Delta_t
       this%phi = this%phi(this%domain%nNodes:1:-1)
       call atx_cr (                          &
            n = this%domain%nNodes            &
            , nz_num = size(worksparse1%A)    &
            , ia = worksparse1%AI             &
            , ja = worksparse1%AJ             &
            , a = worksparse1%A               &
            , x = this%phi                    &
            , w = w                          )
       call atx_cr (                          &
            n = this%domain%nNodes            &
            , nz_num = size(this%MInverse%A)  &
            , ia = this%MInverse%AI           &
            , ja = this%MInverse%AJ           &
            , a = this%MInverse%A             &
            , x = (this%rhs+w)                &
            , w = this%phi                   )
       this%phi = this%phi(this%domain%nNodes:1:-1)
!!$------------------------------------------------------------------
!!$ Boundary conditions
!!$-------------------------------------------------------------------
       Do i = 1, size(this%DirichletNode)
          this%phi(this%DirichletNode(i)) = this%DirichletValue(i)
       End Do
!!$--------------------------------------------------------------------
       call this%computeFlux
       step1 = step1 + 1
    End Do
    Close(Data)
    Close(results)
  End Subroutine TransientMethod

!!$  Subroutine Runge_Kutta4(this)
!!$    Implicit none
!!$    Class(Poisson1DType), Intent(InOut) :: this
!!$  End Subroutine Runge_Kutta4
  
  Subroutine valueJacobianGaussPoints(this)
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    Integer :: nNodes
    Allocate(this%domain%element(iElem)%jacobian(this%nGauss))
    this%domain%element(iElem)%jacobian(:) = 0
    nNodes = this%domain%element(iElem)%nNodes
    Do i = 1, nNodes
       this%domain%element(iElem)%jacobian(:) = this%domain%element(iElem)%jacobian(:) &
            + derivShapeFuncMat(nNodes-1, i, :)*this%domain%x(this%domain%element(iElem)%node(i))
    End Do
  End Subroutine valueJacobianGaussPoints
  
  Subroutine valueJacobianNodePoints(this)
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    Real(dp), Dimension(4) :: u
    Real(dp), Dimension(4) :: dN
    Deallocate(this%domain%element(iElem)%jacobian)
    Allocate(this%domain%element(iElem)%jacobian(this%domain%element(iElem)%nNodes))
    this%domain%element(iElem)%jacobian(:) = 0
    u(1) = -1.d0
    u(4) = 1.d0
    Do i = 1, this%domain%element(iElem)%nNodes
       u(2) = -1.d0 +  2.d0/(this%domain%element(iElem)%nNodes-1)
       u(3) = -1.d0 + 2*2.d0/(this%domain%element(iElem)%nNodes-1)
       dN = this%domain%element(iElem)%dN(u(i),this%domain%element(iElem)%nNodes)
       Do j = 1, this%domain%element(iElem)%nNodes
          this%domain%element(iElem)%jacobian(i) = this%domain%element(iElem)%jacobian(i) &
               + dN(j)*this%domain%x(this%domain%element(iElem)%node(j))
       End Do
    End Do
  End Subroutine valueJacobianNodePoints

  Subroutine computeK(this)
    Implicit none
    Class(poisson1DType), Intent(InOut) :: this
    Integer :: nNodes
    Call this%K%init(nnz = 16*this%domain%nElem, rows = this%domain%nNodes+1)
    Allocate(kVect(this%nGauss))
    Do iElem = 1, this%domain%nElem
       nNodes = this%domain%element(iElem)%nNodes
       Call valueGauss3(this%kCoef &
            , this%domain%element(iElem)%materialIndex &
            , kVect &
            , lower = this%domain%x(this%domain%element(iElem)%node(1)) &
            , upper = this%domain%x(this%domain%element(iElem)%node(nNodes)))
       Do i = 1, nNodes
          Do j = 1, nNodes 
             Kij = gauss(1.d0/this%domain%element(iElem)%jacobian, kVect &
                  , derivShapeFuncMat(nNodes-1, i, :) &
                  , derivShapeFuncMat(nNodes-1, j, :))
             Call this%K%append(value = Kij &
                  , row = this%domain%element(iElem)%node(i) &
                  , col = this%domain%element(iElem)%node(j))
          End Do
       End Do
    End Do
    Call this%K%getSparse
    Deallocate(kVect)
  End Subroutine computeK

 Subroutine computeMLumped(this)
    Implicit none
    Class(poisson1DType), Intent(InOut) :: this
    Integer :: nNodes
    Real(dp) :: factor
    Call this%MLumped%init(nnz = 2*this%domain%nNodes, rows = this%domain%nNodes+1)
    Allocate(rhoVect(this%nGauss),cVect(this%nGauss))
    Do iElem = 1, this%domain%nElem
       nNodes = this%domain%element(iElem)%nNodes
       if (nNodes == 2) then
          factor = 1.5
       else if (nNodes == 3) then
          factor = 2.
       else if (nNodes == 4)
          factor = 2.25
       end if
       Call valueGauss(this%c &
            , this%domain%element(iElem)%materialIndex &
            , cVect &
            , lower = this%domain%x(this%domain%element(iElem)%node(1)) &
            , upper = this%domain%x(this%domain%element(iElem)%node(nNodes)))
       Call valueGauss(this%rho &
            , this%domain%element(iElem)%materialIndex &
            , rhoVect &
            , lower = this%domain%x(this%domain%element(iElem)%node(1)) &
            , upper = this%domain%x(this%domain%element(iElem)%node(nNodes)))
       Do i = 1, nNodes
                Mij = gauss(factor*this%domain%element(iElem)%jacobian, cVect,rhoVect &
                     , shapeFuncMat(nNodes-1, i, :) &
                     , shapeFuncMat(nNodes-1, i, :))
             Call this%MLumped%append(value = Mij &
                  , row = this%domain%element(iElem)%node(i) &
                  , col = this%domain%element(iElem)%node(i))
       End Do
    End Do
    Call this%MLumped%getSparse
    Deallocate(rhoVect, cVect)
  End Subroutine computeMLumped

 Subroutine computeM(this)
    Implicit none
    Class(poisson1DType), Intent(InOut) :: this
    Integer :: nNodes
    Call this%M%init(nnz = 16*this%domain%nElem, rows = this%domain%nNodes+1)
    Allocate(rhoVect(this%nGauss),cVect(this%nGauss))
    Do iElem = 1, this%domain%nElem
       nNodes = this%domain%element(iElem)%nNodes
       Call valueGauss(this%c &
            , this%domain%element(iElem)%materialIndex &
            , cVect &
            , lower = this%domain%x(this%domain%element(iElem)%node(1)) &
            , upper = this%domain%x(this%domain%element(iElem)%node(nNodes)))
       Call valueGauss(this%rho &
            , this%domain%element(iElem)%materialIndex &
            , rhoVect &
            , lower = this%domain%x(this%domain%element(iElem)%node(1)) &
            , upper = this%domain%x(this%domain%element(iElem)%node(nNodes)))
       Do i = 1, nNodes
          Do j = 1, nNodes 
             Mij = gauss(this%domain%element(iElem)%jacobian, cVect,rhoVect &
                  , shapeFuncMat(nNodes-1, i, :) &
                  , shapeFuncMat(nNodes-1, j, :))
             Call this%M%append(value = Mij &
                  , row = this%domain%element(iElem)%node(i) &
                  , col = this%domain%element(iElem)%node(j))
          End Do
        End Do
    End Do
    Call this%M%getSparse
    Deallocate(rhoVect, cVect)
  End Subroutine computeM

  subroutine getInverse(this,n,A,inverse)
    Use GMRES  
    Class(Poisson1DType), Intent(InOut) :: this
    type(sparseType) :: A
    type(sparseType) :: inverse
    real(dp) :: y(n),x(n)
    Integer :: n
    Integer, parameter :: ITR_MAX = 1000
    Real(dp), parameter :: TOL_ABS = 1d-15
    Call inverse%init(nnz = 16*this%domain%nElem, rows = this%domain%nNodes)
    do j = 1,this%domain%nNodes
       y = 0.
       y(j) = 1.
       
       Call pmgmres_ilu_cr(                &
              n = n                        &
            , nz_num = size(A%A)           &
            , ia = A%AI                    &
            , ja = A%AJ                    &
            , a = A%A                      &
            , x = x                        &
            , rhs = y                      &
            , itr_max = ITR_MAX            &
            , mr = n-1                     &
            , tol_abs = TOL_ABS            &
            , tol_rel = 1d-15              &
            , verbose = .false.           )
       do i = 1, n
          if(abs(x(i)).gt.1d-5)then
             Call inverse%append(value = x(i), row = i, col = j)
          end if
       end do   
    end do
    Call inverse%getSparse
  end subroutine getInverse
  
  Subroutine computeRHS(this)
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    Integer :: nNodes
    Allocate(this%rhs(size(this%domain%x)))
    this%rhs = 0
    Allocate(QVect(this%nGauss))
    Do iElem = 1, this%domain%nElem
       nNodes = this%domain%element(iElem)%nNodes
       Call valueGauss3(this%source &
            , this%domain%element(iElem)%sourceIndex &
            , QVect &
            , lower = this%domain%x(this%domain%element(iElem)%node(1)) &
            , upper = this%domain%x(this%domain%element(iElem)%node(nNodes)))
       Do i = 1, this%domain%element(iElem)%nNodes
          fi = gauss(this%domain%element(iElem)%jacobian, QVect &
               , shapeFuncMat(nNodes-1, i, :))
          this%rhs(this%domain%element(iElem)%node(i)) = &
               this%rhs(this%domain%element(iElem)%node(i)) + fi
       End Do
    End Do
  End Subroutine computeRHS

  Subroutine boundaryConditions(this)
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    Integer :: index
    Do i = 1, size(this%DirichletNode)
       index = this%K%AI(this%DirichletNode(i))
       Do while(index < this%K%AI(this%DirichletNode(i)+1)) 
          this%K%A(index) = 0.d0
          If(this%K%AJ(index) .eq. this%DirichletNode(i)) then
             this%K%A(index) = 1.d0
          End If
          index = index + 1
       End Do
       this%rhs(this%DirichletNode(i)) = this%DirichletValue(i)
    End Do
    Do i = 1, size(this%NeumannNode)
       this%rhs(this%NeumannNode(i)) = &
            this%rhs(this%NeumannNode(i)) &
            - this%NeumannValue(i)
    End Do
  End Subroutine boundaryConditions

  Subroutine computeFlux(this)
    use DerivatorMOD
    Implicit none
    Class(Poisson1DType), Intent(InOut) :: this
    real(dp) :: h
    h = this%domain%x(2)-this%domain%x(1) !Temp
    this%flux(1) = getDeriv(this%phi(1:3), h, 'forward')
    Do i = 2, this%domain%nNodes-1
       this%flux(i) = getDeriv(this%phi(i-1:i+1), h, 'centered')
    end Do
    this%flux(this%domain%nNodes) = getDeriv(this%phi(this%domain%nNodes-2:this%domain%nNodes) &
         , h, 'backward')
    this%flux = -1*this%kCoef(1, 1.d0)*this%flux !temp    
  End Subroutine computeFlux
End Module Poisson1DMod


