#ifndef Constitutive_
#define Constitutive_

//2D material model (stress expressions)

#define E4_2D (Ed/(std::pow(Es,4)))
#define E2_2D (Ed/(std::pow(Es,2)))*(2*c-1.0)
#define Eii_2D (Ed/std::pow(Es,2))
#define Eij_2D (Ed/std::pow(Es,2))

//gradient model with nonconvexities (for configurational displacement)
#define PiJ_2D (2*Eii_2D*e1*e1_FiJ + 2*Eij_2D*e6*e6_FiJ + (4*E4_2D*e2*e2*e2 - 4*E2_2D*e2)*e2_FiJ + 2*El*El*Eii_2D*(e2_1*e2_1_FiJ + e2_2*e2_2_FiJ))//2D
#define BetaiJK_2D  (2*El*El*Eii_2D*(e2_1*e2_1_FiJK + e2_2*e2_2_FiJK)) //2D


//3D material model (stress expressions)

#define E4_3D (3*Ed/(2*std::pow(Es,4)))
#define E3_3D (Ed/(std::pow(Es,3)))*(c)
#define E2_3D (3*Ed/(2*std::pow(Es,2)))*(2*c-1.0)
#define Eii_3D (1.5*Ed/std::pow(Es,2))
#define Eij_3D (1.5*Ed/std::pow(Es,2))

//gradient model with nonconvexities (for configurational displacement)
#define PiJ_3D (2*Eii_3D*e1*e1_FiJ + 2*Eij_3D*e4*e4_FiJ + 2*Eij_3D*e5*e5_FiJ + 2*Eij_3D*e6*e6_FiJ + (-2*E2_3D*e2-6*E3_3D*e2*e3+4*E4_3D*e2*(e2*e2+e3*e3))*e2_FiJ + (-2*E2_3D*e3+3*E3_3D*(e3*e3-e2*e2)+4*E4_3D*e3*(e2*e2+e3*e3))*e3_FiJ + 2*El*El*Eii_3D*(e2_1*e2_1_FiJ + e2_2*e2_2_FiJ + e2_3*e2_3_FiJ + e3_1*e3_1_FiJ + e3_2*e3_2_FiJ + e3_3*e3_3_FiJ))
#define BetaiJK_3D  2*El*El*Eii_3D*(e2_1*e2_1_FiJK + e2_2*e2_2_FiJK + e2_3*e2_3_FiJK + e3_1*e3_1_FiJK + e3_2*e3_2_FiJK + e3_3*e3_3_FiJK)

#endif
