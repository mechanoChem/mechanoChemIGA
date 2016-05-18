#define DIM 2
#define GridScale 1.0
#define ADSacado //enables Sacado for AD instead of Adept
#if DIM == 2
	#define numVars 36 //162 in 3D, 36 in 2D for 2*DIM dof
#elif DIM == 3
	#define numVars 162 //162 in 3D, 36 in 2D for 2*DIM dof
#endif

//2D physical parameters
#if DIM == 2
  #define Es 0.1
  #define Ed 1.0
  #define E4 (Ed/(std::pow(Es,4)))
  #define E2 (Ed/(std::pow(Es,2)))*(2*c-1.0)
  #define Eii (Ed/std::pow(Es,2))
  #define Eij (Ed/std::pow(Es,2))
  #define El .025 //**ELambda - constant for gradE.gradE
#elif DIM == 3
  #define Es 0.1
  #define Ed -1.0
  #define E4 (-3*Ed/(2*std::pow(Es,4)))
  #define E3 (-Ed/(std::pow(Es,3)))*(c)
  #define E2 (3*Ed/(2*std::pow(Es,2)))*(2*c-1.0)
  #define Eii (-1.5*Ed/std::pow(Es,2))
  #define Eij (-1.5*Ed/std::pow(Es,2))
  #define El 1 //**ELambda - constant for gradE.gradE
#endif

//material model (stress expressions)
//non-gradient anistopropic St-Venant Kirchoff model (for total displacement)
#define mu 1e5
#define betaC 1e5
#define alphaC 2e5
#if DIM == 2
  #define PiJ ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1])*F[i][J] + 2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J])) //2D
#elif DIM == 3
  #define PiJ ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1]+beta[J][2]*E[2][2])*F[i][J] + 2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J]+F[i][2]*E[2][J])) //3D
#endif
#define BetaiJK (0.0)

//gradient model with nonconvexities (for configurational displacement)
#if DIM == 2
  #define P0iJ (2*Eii*e1*e1_chiiJ + 2*Eij*e6*e6_chiiJ + (4*E4*e2*e2*e2 - 4*E2*e2)*e2_chiiJ + 2*El*El*Eii*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ))//2D
  #define Beta0iJK  (2*El*El*Eii*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK)) //2D
#elif DIM == 3
  #define P0iJ (2*Eii*e1*e1_chiiJ + 2*Eij*e4*e4_chiiJ + 2*Eij*e5*e5_chiiJ + 2*Eij*e6*e6_chiiJ + (2*E2*e2-6*E3*e2*e3+4*E4*e2*(e2*e2+e3*e3))*e2_chiiJ + (2*E2*e3+3*E3*(e3*e3-e2*e2)+4*E4*e3*(e2*e2+e3*e3))*e3_chiiJ + 2*El*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ + e2_3*e2_3_chiiJ + e3_1*e3_1_chiiJ + e3_2*e3_2_chiiJ + e3_3*e3_3_chiiJ))
  #define Beta0iJK  2*El*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK + e2_3*e2_3_chiiJK + e3_1*e3_1_chiiJK + e3_2*e3_2_chiiJK + e3_3*e3_3_chiiJK)
#endif

//boundary conditions
#define uDirichlet 0.1
//other variables
#define NVal 4//**
//time stepping
#define dtVal 0.00001 //** // used to set load parameter..so 0<dtVal<1
#define skipOutput 100
//restart file
#define RESTART_IT 0
#define RESTART_TIME 0

