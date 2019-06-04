Module GidDataMod
  Use tools, Only: dp, sp
  Use Poisson1DMod
  Use Print_All
  
  Implicit none
  Public
  type(Poisson1DType) :: localProblem
  Integer, Parameter :: projectData= 1, project= 2, results= 3
  Integer, Parameter :: Functions= 4
  Integer :: date_time(8)
  Character (len = 250) :: projectName, path
  Real(8), external :: kcoef, source, rho, c

  Interface ReadGid
     Module Procedure ReadGidData
  End Interface ReadGid
  
  Interface WriteGid
     Module Procedure WriteGidData
  End Interface WriteGid

  Contains

  Subroutine ReadGidData(problemOutput)
    Implicit none
    Type(Poisson1DType), Intent(Out) :: problemOutput

    Print*,
    Print*,'Reading project data...'
    Call ReadProjectData()
    Print*,
    Print*,'Reading mesh data...'
    Call InitMesh()
    Print*,
    Print*,'Functions initing...'
    Call InitFunctions()
    Print*,
    Print*,'Loading Boundary Conditions...'
    Call BoundaryConditions()
    problemOutput = localProblem
    Print*, 'READ DATA OK!'
End Subroutine ReadGidData
!------------------------------------------------------------------------
  Subroutine WriteGidData(problemInput)
    Implicit none
    Type(Poisson1DType), Intent(In) :: problemInput
    localProblem = problemInput
    Print*,'Writing results...'
    Select case(localProblem%method)
       Case('Transient_State')
       Call ResetFun()
       Stop
    End Select
    Call GidResults()
    Call ResetFun()
   End Subroutine WriteGidData
!------------------------------------------------------------------------  
  Subroutine ReadProjectData()
  Implicit none
  Open(projectData, file = 'projectData.dat')
  Read(projectData,'(250A)') projectName
  Read(projectData,'(250A)') path
  Close(projectData)
  Print*,
  Print*,'Project Name:', trim(projectName)
  Print*,
  Print*,'Path:', trim(path)
End Subroutine ReadProjectData
!------------------------------------------------------------------------
Subroutine InitMesh()
  Implicit none
  Integer :: i, j, iNode, iElem, Aux(4), nodes

  Open(project, file = trim(projectName)//".dat")
  Do i = 1, 8
  Read(project,*)
  End Do
  Read(project,'(15x,250A)') localProblem%method
  print*,'Analisys Type', localproblem%method
  Select case(localProblem%method)
  Case('Transient_State')
     Read(project,'(8x,250A)') localProblem%SolveMethod
  print*,'Method', localproblem%solvemethod
  Read(project,*)
  Read(project,*) localProblem%t0
  Read(project,*)
  Read(project,*) localProblem%tf
  End Select
  Do i = 1, 4
  Read(project,*)
  End Do
  Read(project,*) localProblem%domain%nElem
  Print*,
  Print*,'Number of Elements:',localProblem%domain%nElem
  Read(project,*)
  Read(project,*) localProblem%domain%nNodes
  Print*,
  Print*,'Number of Nodes:',localProblem%domain%nNodes
  Read(project,*)
  Read(project,*) localProblem%nGauss
  Print*,
  Print*,'Order of Gauss:',localProblem%nGauss
  Do i = 1, 5
  Read(project,*)
  End Do

  Allocate(localProblem%domain%element(localProblem%domain%nElem))
  Aux = 0
  Print*,
  Print*,'Allocated elements, dimension:',size(localProblem%domain%element)
  Print*,
  Print*,'      Element  Material index  source index nodes number     connectivities'
  Do iElem = 1, localProblem%domain%nElem
  Read(project,*) i, localProblem%domain%element(i)%materialIndex, &
                     localProblem%domain%element(i)%sourceIndex, &
                     iNode,(Aux(j),j=1,iNode)
  Select Case(iNode)
     Case(2)
        Call LocalProblem%domain%element(iElem)%init(2)
        localProblem%domain%element(iElem)%node(1:2) = Aux(1:2)
        nodes = 2
     Case(3)
        Call LocalProblem%domain%element(iElem)%init(3)
        localProblem%domain%element(iElem)%node(1) = Aux(1)
        localProblem%domain%element(iElem)%node(2) = Aux(3)
        localProblem%domain%element(iElem)%node(3) = Aux(2)
        nodes = 3
     Case(4)
        Call LocalProblem%domain%element(i)%init(4)
        localProblem%domain%element(i)%node(1:4) = Aux(1:4)
        nodes = 4
  End Select

  Print*, i, localProblem%domain%element(i)%materialIndex, &
             localProblem%domain%element(i)%sourceIndex, &
             nodes,(localProblem%domain%element(iElem)%node(j),j=1,nodes)
  End Do

  Do i = 1, 5
  Read(project,*)
  End Do

  Allocate(localProblem%domain%x(localProblem%domain%nNodes))
  Print*,
  Print*,'Allocated coordinates, dimension:',size(localProblem%domain%x)
  Print*,
  Print*,'Node      coordinate'
  Do iNode = 1, localProblem%domain%nNodes
     Read(project,*) i, localProblem%domain%x(i)
     Print*, i, localProblem%domain%x(i)
  End Do
    
End Subroutine InitMesh
!------------------------------------------------------------------------
Subroutine InitFunctions()
  implicit none
    localProblem%kcoef  => kcoef
    localProblem%source => source
Select case(localProblem%method)
Case('Transient_State')
    localProblem%rho => rho
    localProblem%c   => c
End Select
End Subroutine InitFunctions
!-------------------------------------------------------------------------
Subroutine BoundaryConditions()
  Implicit none
  Integer :: i, nPhi_barra, iPhi_barra, nq_barra, iq_barra, nmat, imat, nsource, isource

  Do i = 1, 4
  Read(project,*)
  End Do
  Read(project,*) nmat
  Read(project,*)
  Read(project,*)
  Do imat = 1, (nmat*3+4)
  read(project,*) 
  End Do
  Do i = 1, 4
  Read(project,*)
  End Do
  Read(project,*) nSource
  Read(project,*)
  Read(project,*)
  Do iSource = 1, nSource
     Read(project,*)
  End Do
  
  Do i = 1, 4
  Read(project,*)
  End Do
  Read(project,*) nPhi_barra
  Read(project,*)
  Read(project,*)

  Allocate(localProblem%DirichletNode(nPhi_barra),localProblem%DirichletValue(nPhi_barra))
  Print*,
  Print*,'Allocated Dirichlet conditions, dimension:',size(localProblem%DirichletNode)
  Print*, 'Node     value'
  Do iPhi_barra = 1, nPhi_barra
  Read(project,*) localProblem%DirichletNode(iphi_barra), localProblem%DirichletValue(iPhi_barra)
  Print*, localProblem%DirichletNode(iphi_barra), localProblem%DirichletValue(iPhi_barra)
  End Do
                          
  Do i = 1, 4
  Read(project,*)
  End Do
  Read(project,*) nq_barra
  Read(project,*)
  Read(project,*)

  Allocate(localProblem%NeumannNode(nq_barra),localProblem%NeumannValue(nq_barra))
  Print*,
  Print*,'Allocated Neumann conditions, dimension:',size(localProblem%NeumannNode)
  Print*, 'Node     value'
  Do iq_barra = 1, nq_barra
  Read(project,*) localProblem%NeumannNode(iq_barra), localProblem%NeumannValue(iq_barra)
  Print*, localProblem%NeumannNode(iq_barra), localProblem%NeumannValue(iq_barra)
  End Do
  Close(project)
End Subroutine BoundaryConditions
!-------------------------------------------------------------------------
Subroutine GidResults()
  Implicit none
  Integer :: iNode
   
  Open(results, file = trim(projectName)//'.flavia.res')
  write(results,*)      'GiD Post Result File 1.0'
  write(results,*)      'Result "Temperature"',' "',trim(projectName),'" ','1 Scalar OnNodes'
  write(results,*)      'Values' 
  Do iNode = 1, localProblem%domain%nNodes
  Write(results,*)      iNode, localProblem%phi(iNode)
  End Do
  write(results,*)      'End Values'
  Write(results,*)
  write(results,*)      'Result "Flux"',' "',trim(projectName),'" ','1 Scalar OnNodes'
  write(results,*)      'Values'
  Do iNode = 1, localProblem%domain%nNodes
  Write(results,*)      iNode, localProblem%flux(iNode)
  End Do
  write(results,*)      'End Values'
  Close(results)
End Subroutine GidResults
!----------------------------------------------------------------------------
Subroutine ResetFun()
  Implicit none
  Open(Functions, file = trim(path)//'/Poisson1D/Source/functions.f90')
  write(Functions,*)    'Function kCoef(i,Xvalue)' 
  write(Functions,*)    'Implicit none'
  write(Functions,*)    'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)    'Real(8) :: kcoef'
  write(Functions,*)    'Integer, Intent(In) :: i'
  write(Functions,*)    'end function kCoef'
  write(Functions,*)
  write(Functions,*)
  write(Functions,*)   'Function Source(i,Xvalue)' 
  write(Functions,*)   'Implicit none'
  write(Functions,*)   'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)   'Real(8) :: Source'
  write(Functions,*)   'Integer, Intent(In) :: i'
  write(Functions,*)   'End Function Source'
  write(Functions,*) 
  write(Functions,*)
  write(Functions,*)    'Function rho(i,Xvalue)' 
  write(Functions,*)    'Implicit none'
  write(Functions,*)    'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)    'Real(8) :: rho'
  write(Functions,*)    'Integer, Intent(In) :: i'
  write(Functions,*)    'end function rho'
  write(Functions,*)
  write(Functions,*)
  write(Functions,*)    'Function c(i,Xvalue)' 
  write(Functions,*)    'Implicit none'
  write(Functions,*)    'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)    'Real(8) :: c'
  write(Functions,*)    'Integer, Intent(In) :: i'
  write(Functions,*)    'end function c'
  write(Functions,*)
  close(Functions)

  Call date_and_time(VALUES=date_time)
  Print*,'Finish Poisson 1D'
  Print*,'Date:',date_time(3),"/",date_time(2),"/", date_time(1)
  Print*,'Hour:',date_time(5),":",date_time(6),":",date_time(7)
End Subroutine ResetFun
!----------------------------------------------------------------------------
End Module GidDataMod
