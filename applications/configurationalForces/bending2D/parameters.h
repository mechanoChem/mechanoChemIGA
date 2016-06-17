/*#define DIM 2
#define GridScale 1.0
#define ADSacado //enables Sacado for AD instead of Adept
#define numVars 36 //162 in 3D, 36 in 2D for 2*DIM dof
*/
//2D physical parameters
#define Es 0.1
#define Ed 1.0
#define E4 (Ed/(std::pow(Es,4)))
#define E2 (Ed/(std::pow(Es,2)))*(2*c-1.0)
#define Eii (Ed/std::pow(Es,2))
#define Eij (Ed/std::pow(Es,2))
#define El .025 //**ELambda - constant for gradE.gradE

//material model (stress expressions)
//non-gradient anistopropic St-Venant Kirchoff model (for total displacement)
#define mu 1e5
#define betaC 1e5
#define alphaC 2e5
#define PiJ ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1])*F[i][J] + 2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J])) //2D
#define BetaiJK (0.0)

//gradient model with nonconvexities (for configurational displacement)
#define P0iJ (2*Eii*e1*e1_chiiJ + 2*Eij*e6*e6_chiiJ + (4*E4*e2*e2*e2 - 4*E2*e2)*e2_chiiJ + 2*El*El*Eii*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ))//2D
#define Beta0iJK  (2*El*El*Eii*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK)) //2D

//anisotropy
#define alphaI (alphaC*Lambda[I]) //Stiffens with elongation
#define d_alphaL (alphaC/Lambda[L])
/*
//boundary conditions
#define uDirichlet 0.1
//other variables
#define NVal 5//**
//time stepping
#define dtVal 0.00001 //** // used to set load parameter..so 0<dtVal<1
#define skipOutput 1
//restart file
#define RESTART_IT 0
#define RESTART_TIME 0
*/
