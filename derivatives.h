#include <math.h>

#define DIM 3

inline PetscReal EbyF(unsigned int A, unsigned int B, unsigned int i, unsigned int J){
  return 0.5*(F[i][A]*((PetscReal)(B==J))+F[i][B]*((PetscReal)(A==J)));
}

inline PetscReal dEbyF(unsigned int A, unsigned int B, unsigned int C, unsigned int i, unsigned int J){
  return 0.5*(dF[i][A][C]*((PetscReal)(B==J))+dF[i][B][C]*((PetscReal)(A==J)));
}

inline PetscReal dEbydF(unsigned int A, unsigned int B, unsigned int C, unsigned int i, unsigned int J, unsigned int K){
  return 0.5*(F[i][A]*((PetscReal)(B==J))+F[i][B]*((PetscReal)(A==J)))*((double) (C==K));
}

inline PetscReal EbyFF(unsigned int A, unsigned int B, unsigned int i, unsigned int J, unsigned int l, unsigned int M){
  return 0.5*((double) (i==L)*((A==M)*(B==J) + (B==M)*(A==J))); 
}

inline PetscReal dEbyFdF(unsigned int A, unsigned int B, unsigned int C, unsigned int i, unsigned int J, unsigned int K, unsigned int l, unsigned int M){
  return 0.5*((double) (i==L)*(C==K)*((A==M)*(B==J) + (B==M)*(A==J))); 
}

void ebyF(PetscReal* arraye_F){
  PetscReal (*e_F)[DIM][DIM]=(PetscReal (*)[DIM][DIM]) arraye_F;
  for (unsigned int i=0; i<DIM; i++)
    for (unsigned int J=0; J<DIM; J++){
      e_F[0][i][J]=(1.0/std::sqrt(3.0))*(EbyF(0,0,i,J)+EbyF(1,1,i,J)+EbyF(2,2,i,J));   //alpha=1
      e_F[1][i][J]=(1.0/std::sqrt(2.0))*(EbyF(0,0,i,J)-EbyF(1,1,i,J)); //alpha=2
      e_F[2][i][J]=(1.0/std::sqrt(6.0))*(EbyF(0,0,i,J)+EbyF(1,1,i,J)-2*EbyF(2,2,i,J));   //alpha=3
      e_F[3][i][J]=EbyF(1,2,i,J); //alpha=4 
      e_F[4][i][J]=EbyF(2,0,i,J); //alpha=5 
      e_F[5][i][J]=EbyF(0,1,i,J); //alpha=6 
    }
}

void debyF(PetscReal* arrayde_F){
  PetscReal (*de_F)[DIM][DIM][DIM]=(PetscReal (*)[DIM][DIM][DIM]) arrayde_F;
  for (unsigned int b=0; b<DIM; b++)
    for (unsigned int i=0; i<DIM; i++)
      for (unsigned int J=0; J<DIM; J++){
	de_F[0][b][i][J]=(1.0/std::sqrt(3.0))*(dEbyF(0,0,b,i,J)+dEbyF(1,1,b,i,J)+dEbyF(2,2,b,i,J));   //alpha=1
	de_F[1][b][i][J]=(1.0/std::sqrt(2.0))*(dEbyF(0,0,b,i,J)-dEbyF(1,1,b,i,J)); //alpha=2
	de_F[2][b][i][J]=(1.0/std::sqrt(6.0))*(dEbyF(0,0,b,i,J)+dEbyF(1,1,b,i,J)-2*dEbyF(2,2,b,i,J));   //alpha=3
	de_F[3][b][i][J]=dEbyF(1,2,b,i,J); //alpha=4 
	de_F[4][b][i][J]=dEbyF(2,0,b,i,J); //alpha=5 
	de_F[5][b][i][J]=dEbyF(0,1,b,i,J); //alpha=6 
      }
}

void debydF(PetscReal* arrayde_dF){
  PetscReal (*de_dF)[DIM][DIM][DIM][DIM]=(PetscReal (*)[DIM][DIM][DIM][DIM]) arrayde_dF;
  for (unsigned int b=0; b<DIM; b++)
    for (unsigned int i=0; i<DIM; i++)
      for (unsigned int J=0; J<DIM; J++)
	for (unsigned int K=0; K<DIM; K++){
	  de_dF[0][b][i][J][K]=(1.0/std::sqrt(3.0))*(dEbydF(0,0,b,i,J,K)+dEbydF(1,1,b,i,J,K)+dEbydF(2,2,b,i,J,K));   //alpha=1
	  de_dF[1][b][i][J][K]=(1.0/std::sqrt(2.0))*(dEbydF(0,0,b,i,J,K)-dEbydF(1,1,b,i,J,K)); //alpha=2
	  de_dF[2][b][i][J][K]=(1.0/std::sqrt(6.0))*(dEbydF(0,0,b,i,J,K)+dEbydF(1,1,b,i,J,K)-2*dEbydF(2,2,b,i,J,K));   //alpha=3
	  de_dF[3][b][i][J][K]=dEbydF(1,2,b,i,J,K); //alpha=4 
	  de_dF[4][b][i][J][K]=dEbydF(2,0,b,i,J,K); //alpha=5 
	  de_dF[5][b][i][J][K]=dEbydF(0,1,b,i,J,K); //alpha=6 
	}
}

void ebyFF(PetscReal* arraye_FF){
  PetscReal (*e_FF)[DIM][DIM][DIM][DIM]=(PetscReal (*)[DIM][DIM][DIM][DIM]) arraye_FF;
  for (unsigned int i=0; i<DIM; i++)
    for (unsigned int J=0; J<DIM; J++)
      for (unsigned int l=0; l<DIM; l++)
	for (unsigned int M=0; M<DIM; M++){
	  e_FF[0][i][J][l][M]=(1.0/std::sqrt(3.0))*(EbyFF(0,0,i,J,l,M)+EbyFF(1,1,i,J,l,M)+EbyFF(2,2,i,J,l,M));   //alpha=1
	  e_FF[1][i][J][l][M]=(1.0/std::sqrt(2.0))*(EbyFF(0,0,i,J,l,M)-EbyFF(1,1,i,J,l,M)); //alpha=2
	  e_FF[2][i][J][l][M]=(1.0/std::sqrt(6.0))*(EbyFF(0,0,i,J,l,M)+EbyFF(1,1,i,J,l,M)-2*EbyFF(2,2,i,J,l,M));   //alpha=3
	  e_FF[3][i][J][l][M]=EbyFF(1,2,i,J,l,M); //alpha=4 
	  e_FF[4][i][J][l][M]=EbyFF(2,0,i,J,l,M); //alpha=5 
	  e_FF[5][i][J][l][M]=EbyFF(0,1,i,J,l,M); //alpha=6 
	}
}

void debyFdF(PetscReal* arrayde_FdF){
  PetscReal (*de_FdF)[DIM][DIM][DIM][DIM][DIM][DIM]=(PetscReal (*)[DIM][DIM][DIM][DIM][DIM][DIM]) arrayde_FdF;
  for (unsigned int b=0; b<DIM; b++)
    for (unsigned int i=0; i<DIM; i++)
      for (unsigned int J=0; J<DIM; J++)
	for (unsigned int K=0; K<DIM; K++)
	  for (unsigned int l=0; l<DIM; l++)
	    for (unsigned int M=0; M<DIM; M++){
	      de_FdF[0][b][i][J][K][l][M]=(1.0/std::sqrt(3.0))*(dEbyFdF(0,0,b,i,J,K,l,M)+dEbyFdF(1,1,b,i,J,K,l,M)+dEbyFdF(2,2,b,i,J,K,l,M));   //alpha=1
	      de_FdF[1][b][i][J][K][l][M]=(1.0/std::sqrt(2.0))*(dEbyFdF(0,0,b,i,J,K,l,M)-dEbyFdF(1,1,b,i,J,K,l,M)); //alpha=2
	      de_FdF[2][b][i][J][K][l][M]=(1.0/std::sqrt(6.0))*(dEbyFdF(0,0,b,i,J,K,l,M)+dEbyFdF(1,1,b,i,J,K,l,M)-2*dEbyFdF(2,2,b,i,J,K,l,M));   //alpha=3
	      de_FdF[3][b][i][J][K][l][M]=dEbyFdF(1,2,b,i,J,K,l,M); //alpha=4 
	      de_FdF[4][b][i][J][K][l][M]=dEbyFdF(2,0,b,i,J,K,l,M); //alpha=5 
	      de_FdF[5][b][i][J][K][l][M]=dEbyFdF(0,1,b,i,J,K,l,M); //alpha=6 
	    }
}

void phibye(PetscReal* arrayphi_e, PetscReal* arrayE){
  PetscReal (*phi_e)=(PetscReal (*)) arrayphi_e;
  PetscReal (*E)[DIM]=(PetscReal (*)[DIM]) arrayE;
  PetscReal e1=(1.0/std::sqrt(3.0))*(E[0][0]+E[1][1]+E[2][2]);
  PetscReal e2=(1.0/std::sqrt(2.0))*(E[0][0]-E[1][1]);
  PetscReal e3=(1.0/std::sqrt(6.0))*(E[0][0]+E[1][1]-2*E[2][2]);
  PetscReal e4=E[1][2];
  PetscReal e5=E[2][0];
  PetscReal e6=E[0][1];
  
  phi_e[0]=2*A*e1; //alpha=1
  phi_e[1]=2*B*e2-6*C*e2*e3+4*D*(e2*e2+e3*e3)*e2; //alpha=2
  phi_e[2]=2*B*e3+3*C*(e3*e3-e2*e2)+4*D*(e2*e2+e3*e3)*e3; //alpha=3
  phi_e[3]=2*E*e4; //alpha=4
  phi_e[4]=2*E*e5; //alpha=5
  phi_e[5]=2*E*e6; //alpha=6
}

void phibyee(PetscReal* arrayphi_ee, PetscReal* arrayE){
  PetscReal (*phi_ee) [6]=(PetscReal (*) [6]) arrayphi_ee;
  PetscReal (*E)[DIM]=(PetscReal (*)[DIM]) arrayE;
  PetscReal e1=(1.0/std::sqrt(3.0))*(E[0][0]+E[1][1]+E[2][2]);
  PetscReal e2=(1.0/std::sqrt(2.0))*(E[0][0]-E[1][1]);
  PetscReal e3=(1.0/std::sqrt(6.0))*(E[0][0]+E[1][1]-2*E[2][2]);
  PetscReal e4=E[1][2];
  PetscReal e5=E[2][0];
  PetscReal e6=E[0][1];
  
  phi_ee[0]=2*A*e1; //alpha=1
  phi_ee[1]=2*B*e2-6*C*e2*e3+4*D*(e2*e2+e3*e3)*e2; //alpha=2
  phi_ee[2]=2*B*e3+3*C*(e3*e3-e2*e2)+4*D*(e2*e2+e3*e3)*e3; //alpha=3
  phi_ee[3]=2*E*e4; //alpha=4
  phi_ee[4]=2*E*e5; //alpha=5
  phi_ee[5]=2*E*e6; //alpha=6

  for (unsigned int a=0; a<6; a++)
    for (unsigned int b=0; b<6; b++){
      if ((a==0) && (b==0)) phi_ee[a][b]=2*A;
      else if ((a==1) && (b==1)) phi_ee[a][b]=2*B-6*C*e3+12*D*e2*e2;
      else if ((a==2) && (b==2)) phi_ee[a][b]=2*B+6*C*e3+12*D*e3*e3;
      else if ((a==1) && (b==2)) phi_ee[a][b]=-6*C*e2+8*D*e3*e2;
      else if ((a==2) && (b==1)) phi_ee[a][b]=-6*C*e2+8*D*e3*e2;
      else if ((a==3) && (b==3)) phi_ee[a][b]=2*E;
      else if ((a==4) && (b==4)) phi_ee[a][b]=2*E;
      else if ((a==5) && (b==5)) phi_ee[a][b]=2*E; 
      else phi_ee[a][b]=0.0;
    }
}

void psibyde(PetscReal* arraypsi_de, PetscReal* arraydE){
  PetscReal (*psi_de)[DIM]=(PetscReal (*) [DIM]) arraypsi_de;
  PetscReal (*dE)[DIM][DIM]=(PetscReal (*)[DIM][DIM]) arraydE;

  PetscReal e21=(1.0/std::sqrt(2.0))*(dE[0][0][0]-dE[1][1][0]);
  PetscReal e22=(1.0/std::sqrt(2.0))*(dE[0][0][1]-dE[1][1][1]);
  PetscReal e23=(1.0/std::sqrt(2.0))*(dE[0][0][2]-dE[1][1][2]);
  PetscReal e31=(1.0/std::sqrt(6.0))*(dE[0][0][0]+dE[1][1][0]-2*dE[2][2][0]);
  PetscReal e32=(1.0/std::sqrt(6.0))*(dE[0][0][1]+dE[1][1][1]-2*dE[2][2][1]);
  PetscReal e33=(1.0/std::sqrt(6.0))*(dE[0][0][2]+dE[1][1][2]-2*dE[2][2][2]);
  
  psi_de[0][0]=2*g*e21; //alpha=2, beta=1
  psi_de[0][1]=2*g*e22; //alpha=2, beta=2
  phi_de[0][2]=2*g*e23; //alpha=2, beta=3
  psi_de[1][0]=2*g*e31; //alpha=3, beta=1
  psi_de[1][1]=2*g*e32; //alpha=3, beta=2
  phi_de[1][2]=2*g*e33; //alpha=3, beta=3
}

void psibydede(PetscReal* arraypsi_dede){
  PetscReal (*psi_dede)[DIM][2][DIM]=(PetscReal (*) [DIM][2][DIM]) arraypsi_dede;

  for (unsigned int a=0; a<2; a++)
    for (unsigned int b=0; b<DIM; b++)
      for (unsigned int g=0; g<2; g++)
	for (unsigned int d=0; d<DIM; d++){
	  if ((a==g) && (b==d)) psi_dede[a][b][g][d]=2*g;
	  else  psi_dede[a][b][g][d]=0.0;
	}
}
