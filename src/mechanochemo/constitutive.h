#ifndef Constitutive_
#define Constitutive_

//2D material model

#define C4_2D (-16*Cd/std::pow(Cs,4))
#define C3_2D (32*Cd/std::pow(Cs,3))
#define C2_2D (-16*Cd/std::pow(Cs,2))
#define E4_2D (-Ed/std::pow(Es,4))
#define E3_2D 0.0
#define E2_2D (2*Ed/std::pow(Es,2))*((2*c-Cs)/Cs)
#define E4_c_2D 0.0
#define E3_c_2D 0.0
#define E2_c_2D (2*Ed/std::pow(Es,2))*(2.0/Cs)
#define E4_cc_2D 0.0
#define E3_cc_2D 0.0
#define E2_cc_2D 0.0
#define Eii_2D (-2*Ed/std::pow(Es,2))
#define Eij_2D (-2*Ed/std::pow(Es,2))

//stress
#define PiJ_2D (2*Eii_2D*e1*e1_FiJ + 2*Eij_2D*e6*e6_FiJ + (4*E4_2D*e2*e2*e2 + 3*E3_2D*e2*e2 + 2*E2_2D*e2)*e2_FiJ + (El*e2_1+Gl*cx[0]/2)*e2_1_FiJ + (El*e2_2+Gl*cx[1]/2)*e2_2_FiJ)
#define BetaiJK_2D ((El*e2_1+Gl*cx[0]/2)*e2_1_FiJK + (El*e2_2+Gl*cx[1]/2)*e2_2_FiJK)

//chemistry
#define mu_c_2D (12*C4_2D*c*c+6*C3_2D*c+2*C2_2D+E4_cc_2D*e2*e2*e2*e2+E3_cc_2D*e2*e2*e2+E2_cc_2D*e2*e2)
#define mu_e2_2D (4*E4_c_2D*e2*e2*e2+3*E3_c_2D*e2*e2+2*E2_c_2D*e2)

//3D material model

#define C4_3D (-16*Cd/std::pow(Cs,4))
#define C3_3D (32*Cd/std::pow(Cs,3))
#define C2_3D (-16*Cd/std::pow(Cs,2))
#define E4_3D (-3*Ed/(2*std::pow(Es,4)))
#define E3_3D (Ed/(std::pow(Es,3)))*(c/Cs)
#define E2_3D (3*Ed/(2*std::pow(Es,2)))*(2*c/Cs-1.0)
#define E4_c_3D 0.0
#define E3_c_3D (Ed/(std::pow(Es,3)))*(1.0/Cs)
#define E2_c_3D (3*Ed/(2*std::pow(Es,2)))*(2.0/Cs)
#define E4_cc_3D 0.0
#define E3_cc_3D 0.0
#define E2_cc_3D 0.0
#define Eii_3D (-1.5*Ed/std::pow(Es,2))
#define Eij_3D (-1.5*Ed/std::pow(Es,2))

//stress
#define PiJ_3D (2*Eii_3D*e1*e1_FiJ + 2*Eij_3D*e4*e4_FiJ + 2*Eij_3D*e5*e5_FiJ + 2*Eij_3D*e6*e6_FiJ + (2*E2_3D*e2-6*E3_3D*e2*e3+4*E4_3D*e2*(e2*e2+e3*e3))*e2_FiJ + (2*E2_3D*e3+3*E3_3D*(e3*e3-e2*e2)+4*E4_3D*e3*(e2*e2+e3*e3))*e3_FiJ + 2*El*(e2_1*e2_1_FiJ + e2_2*e2_2_FiJ + e2_3*e2_3_FiJ + e3_1*e3_1_FiJ + e3_2*e3_2_FiJ + e3_3*e3_3_FiJ))
#define BetaiJK_3D  2*El*(e2_1*e2_1_FiJK + e2_2*e2_2_FiJK + e2_3*e2_3_FiJK + e3_1*e3_1_FiJK + e3_2*e3_2_FiJK + e3_3*e3_3_FiJK)
//chemistry
#define mu_c_3D (12*C4_3D*c*c+6*C3_3D*c+2*C2_3D) // E4_cc, E3_cc, E2_cc are zero so skipping those terms 
#define mu_e2_3D (4*E4_c_3D*e2*(e2*e2+e3*e3) - 6*E3_c_3D*e3*e2 + 2*E2_c_3D*e2)
#define mu_e3_3D (4*E4_c_3D*e3*(e2*e2+e3*e3) + 3*E3_c_3D*(e3*e3-e2*e2) + 2*E2_c_3D*e3)

#endif
