#ifndef Expressions_
#define Expressions_

//3D physical parameters
#define Es 0.1
#define Ed -1.0
#define E4 (-3*Ed/(2*std::pow(Es,4)))
#define E3 (-Ed/(std::pow(Es,3)))*(c)
#define E2 (3*Ed/(2*std::pow(Es,2)))*(2*c-1.0)
#define Eii (-1.5*Ed/std::pow(Es,2))
#define Eij (-1.5*Ed/std::pow(Es,2))
#define El 1 //**ELambda - constant for gradE.gradE

//material model (stress expressions)
//non-gradient anistopropic St-Venant Kirchoff model (for total displacement)
#define mu 1e5
#define betaC 1e5
#define alphaC 2e5
#define PiJ ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1]+beta[J][2]*E[2][2])*F[i][J] + 2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J]+F[i][2]*E[2][J])) //3D
#define BetaiJK (0.0)

//gradient model with nonconvexities (for configurational displacement)
#define P0iJ (2*Eii*e1*e1_chiiJ + 2*Eij*e4*e4_chiiJ + 2*Eij*e5*e5_chiiJ + 2*Eij*e6*e6_chiiJ + (2*E2*e2-6*E3*e2*e3+4*E4*e2*(e2*e2+e3*e3))*e2_chiiJ + (2*E2*e3+3*E3*(e3*e3-e2*e2)+4*E4*e3*(e2*e2+e3*e3))*e3_chiiJ + 2*El*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ + e2_3*e2_3_chiiJ + e3_1*e3_1_chiiJ + e3_2*e3_2_chiiJ + e3_3*e3_3_chiiJ))
#define Beta0iJK  2*El*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK + e2_3*e2_3_chiiJK + e3_1*e3_1_chiiJK + e3_2*e3_2_chiiJK + e3_3*e3_3_chiiJK)

#endif
