#ifndef Expressions_
#define Expressions_

//2D material model (stress expressions)

#define E4_2D (user->Ed/(std::pow(user->Es,4)))
#define E2_2D (user->Ed/(std::pow(user->Es,2)))*(2*c-1.0)
#define Eii_2D (user->Ed/std::pow(user->Es,2))
#define Eij_2D (user->Ed/std::pow(user->Es,2))

//non-gradient anistopropic St-Venant Kirchoff model (for total displacement)
#define PiJ_2D ((alpha[J]-2*user->mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1])*F[i][J] + 2*user->mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J])) //2D
#define BetaiJK_2D (0.0)

//gradient model with nonconvexities (for configurational displacement)
#define P0iJ_2D (2*Eii_2D*e1*e1_chiiJ + 2*Eij_2D*e6*e6_chiiJ + (4*E4_2D*e2*e2*e2 - 4*E2_2D*e2)*e2_chiiJ + 2*user->El*user->El*Eii_2D*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ))//2D
#define Beta0iJK_2D  (2*user->El*user->El*Eii_2D*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK)) //2D


//3D material model (stress expressions)

#define E4_3D (3*user->Ed/(2*std::pow(user->Es,4)))
#define E3_3D (user->Ed/(std::pow(user->Es,3)))*(c)
#define E2_3D (3*user->Ed/(2*std::pow(user->Es,2)))*(2*c-1.0)
#define Eii_3D (1.5*user->Ed/std::pow(user->Es,2))
#define Eij_3D (1.5*user->Ed/std::pow(user->Es,2))

//non-gradient anistopropic St-Venant Kirchoff model (for total displacement)
#define PiJ_3D ((alpha[J]-2*user->mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1]+beta[J][2]*E[2][2])*F[i][J] + 2*user->mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J]+F[i][2]*E[2][J])) //3D
#define BetaiJK_3D (0.0)

//gradient model with nonconvexities (for configurational displacement)
#define P0iJ_3D (2*Eii_3D*e1*e1_chiiJ + 2*Eij_3D*e4*e4_chiiJ + 2*Eij_3D*e5*e5_chiiJ + 2*Eij_3D*e6*e6_chiiJ + (-2*E2_3D*e2-6*E3_3D*e2*e3+4*E4_3D*e2*(e2*e2+e3*e3))*e2_chiiJ + (-2*E2_3D*e3+3*E3_3D*(e3*e3-e2*e2)+4*E4_3D*e3*(e2*e2+e3*e3))*e3_chiiJ + 2*user->El*user->El*Eii_3D*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ + e2_3*e2_3_chiiJ + e3_1*e3_1_chiiJ + e3_2*e3_2_chiiJ + e3_3*e3_3_chiiJ))
#define Beta0iJK_3D  2*user->El*user->El*Eii_3D*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK + e2_3*e2_3_chiiJK + e3_1*e3_1_chiiJK + e3_2*e3_2_chiiJK + e3_3*e3_3_chiiJK)

#endif
