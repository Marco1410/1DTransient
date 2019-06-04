###########################################################
                PROGRAM POISSON'S ANALYSIS
###########################################################

##################### PROBLEM DATA ########################

Problem Name: *gendata(Problem_Name)
Solver Type: *gendata(Solver_Type)
Solver Method: *gendata(Solver_Method)
*if(strcmp(gendata(Solver_Method),"Transient_State")==0)
Method: *\
*if(strcmp(gendata(Explicit),"YES")==0)
Explicit
Initial Time: 
*gendata(Initial_Time)
Final Time: 
*gendata(Final_Time)
*endif
*if(strcmp(gendata(Crank-Nicolson),"YES")==0)
Crank-Nicolson
Initial Time: 
*gendata(Initial_Time)
Final Time: 
*gendata(Final_Time)
*endif
*if(strcmp(gendata(Galerkin),"YES")==0)
Galerkin
Initial Time: 
*gendata(Initial_Time)
Final Time: 
*gendata(Final_Time)
*endif
*if(strcmp(gendata(Implicit),"YES")==0)
Implicit
Initial Time: 
*gendata(Initial_Time)
Final Time: 
*gendata(Final_Time)
*endif
*if(strcmp(gendata(Runge_Kutta4),"YES")==0)
Runge_Kutta4
Initial Time: 
*gendata(Initial_Time)
Final Time: 
*gendata(Final_Time)
*endif
*endif

########################## Mesh ###########################

Elements Number: 
*nelem
Nodes Number: 
*npoin
Gauss Order:
*GenData(Gauss_Order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Element List:
      Element   Material    Source      Nodes        Conectivities
*Set Cond Source *elems
*loop elems
*format "%10i%10i%10i%10i%10i%10i"
*elemsnum *elemsmat *cond(Source_Number) *ElemsNnode *elemsconec
*end elems

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Coordinates:
  Node        X               
*set elems(all)
*loop nodes
*format "%5i%14.5e%14.5e"
*NodesNum *NodesCoord(1,real)
*end nodes

######################## Materials ########################

Materials Number:
*nmats
Materials List:
Material     Thermal conductivity
*loop materials
*matnum *matprop(Thermal_Conductivity)
*end

Material     Density
*loop materials
*matnum *matprop(Density)
*end

Material     Heat specific
*loop materials
*matnum *matprop(specific_Heat)
*end

######################### Sources #########################

*Set Cond Source *elems
Sources Number: 
*Gendata(Source_Number,int)
Conditions List:
      Source    Function
*for(i=1;i<=Gendata(Source_Number,int);i=i+1))
*loop elems
*if(i==cond(Source_Number,int)) 
*cond(Source_Number) *cond(Source) 
*break
*endif
*end
*end for

####################### Temperatures ######################

*Set Cond Fix_Temperature *nodes
Conditions Number: 
*condnumentities
Conditions List: 
       Node        Temperature
*loop nodes *OnlyInCond
*format "%10i%14i"
*NodesNum *cond(Temperature)
*end

########################## Fluxs ##########################

*Set Cond Fix_Flux *nodes
Conditions Number: 
*condnumentities
Conditions List:
       Node     Flux
*loop nodes *OnlyInCond
*format "%10i%8i"
*NodesNum *cond(Flux)
*end

###########################################################