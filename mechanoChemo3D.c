#include <math.h> 
extern "C" {
#include "petiga.h"
}
#include "fields.h"
#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> doubleAD;

#define DIM 2
#define EXPLICIT 1
#define PI 3.14159265

typedef struct {
  IGA iga;
  PetscBool IMPLICIT;
  PetscReal Cs, Cd, C4, C2, CLambda;
  PetscReal Es, Ed, E4, E3, E2, Eii, Eij, Eg, Eh, El;
  PetscReal D, flux;
  PetscReal dt;
  PetscReal gamma, C, he;
  PetscReal cbar;
} AppCtx;

#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
				PetscReal shift,const PetscScalar *V,
				PetscReal t,const T *U,
				PetscReal t0,const PetscScalar *U0,
				T *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;
  PetscBool IMPLICIT = user->IMPLICIT;

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;
  
  //concentration field variable
  T c, cx[DIM], cxx[DIM][DIM]; PetscReal c0;
  computeField<T,DIM,DIM+1>(SCALAR,0,p,U,&c,&cx[0],&cxx[0][0]);
  computeField<PetscReal,DIM,DIM+1>(SCALAR,0,p,U0,&c0);

  //displacement field variables
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DIM+1>(VECTOR,1,p,U,&u[0],&ux[0][0],&uxx[0][0][0]);

  //problem parameters
  PetscReal Cs=user->Cs;
  PetscReal Cd=user->Cd;
  PetscReal C4=user->C4;
  PetscReal C2=user->C2;
  PetscReal CLambda=user->CLambda;
  PetscReal Es=user->Es;
  PetscReal Ed=user->Ed;
  PetscReal E4=user->E4;
  PetscReal E3=user->E3;
  PetscReal E2=user->E2;
  PetscReal Eii=user->Eii;
  PetscReal Eij=user->Eij;
  PetscReal Eg=user->Eg;
  PetscReal Eh=user->Eh;
  PetscReal El=user->El;
  PetscReal D=user->D;
  PetscReal flux=user->flux;
  PetscReal gamma=user->gamma;
  PetscReal C=user->C;
  PetscReal he=user->he;
  PetscReal dt=user->dt;
  
  //Compute F
  T F[3][3], dF[3][3][3];
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      if ((i<DIM) && (j<DIM))  F[i][j]=(i==j)+ux[i][j];
      else  F[i][j]=(i==j);
      for (unsigned int k=0; k<3; k++) {
	if (((i<DIM) && (j<DIM)) && (k<DIM)) dF[i][j][k]=uxx[i][j][k];
	else dF[i][j][k]=0.0;
      }
    }
  }
  //Compute strain metric, E  (E=0.5*(F^T*F-I))
  T E[3][3];
  for (unsigned int i=0; i<3; i++){
    for (unsigned int j=0; j<3; j++){
      E[i][j] = -0.5*(i==j);
      for (unsigned int k=0; k<3; k++){
	E[i][j] += 0.5*F[k][i]*F[k][j];
      }
    }
  }
  //new strain metrics
  T e1=(E[0][0]+E[1][1]+E[2][2])/sqrt(3.0);
  T e2=(E[0][0]-E[1][1])/sqrt(2.0);
  T e3=(E[0][0]+E[1][1]-2*E[2][2])/sqrt(6.0);
  T e4=E[1][2], e5=E[2][0], e6=E[0][1];

  T e2_1=0.0, e2_2=0.0, e2_3=0.0, e3_1=0.0, e3_2=0.0, e3_3=0.0;
  for (unsigned int i=0; i<DIM; ++i){
    e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0])/sqrt(2.0);
    e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1])/sqrt(2.0);
    e2_3+=(F[i][0]*dF[i][0][2]-F[i][1]*dF[i][1][2])/sqrt(2.0);
    e3_1+=(F[i][0]*dF[i][0][0]+F[i][1]*dF[i][1][0]-2*F[i][2]*dF[i][2][0])/sqrt(6.0);
    e3_2+=(F[i][0]*dF[i][0][1]+F[i][1]*dF[i][1][1]-2*F[i][2]*dF[i][2][1])/sqrt(6.0);
    e3_3+=(F[i][0]*dF[i][0][2]+F[i][1]*dF[i][1][2]-2*F[i][2]*dF[i][2][2])/sqrt(6.0);
  }
  
  //phi=E4(e2^2+e3^2)^2+E3*e3*(e3^2-3*e2^2)+E2*(e2^2+e3^2)+Eii*(e1^2)+Eij*(e4^2+e5^2+e6^2)
  //   +Eg(e2_1^2+e2_2^2+(e3_1^2+e3_2^2)/3+ (2/sqrt(3))*(e2_1*e3_1-e2_2*e3_2))
  //   +Eh(e3_1^2+e3_2^2-sqrt(3)*(e2_1*e3_1-e2_2*e3_2))
  //compute P and Beta
  T P[DIM][DIM], Beta[DIM][DIM][DIM];
  T minusC=(Cs-c)/(2.0*Cs), plusC=(Cs+c)/(2.0*Cs);
  T E2c=E2*minusC,  E3c=E3*plusC, E4c=E4*plusC;
  PetscPrintf(PETSC_COMM_WORLD,"\nminusC: %6.3e, plusC: %6.3e",minusC.val(),plusC.val());
  for (unsigned int i=0; i<DIM; ++i){
    for (unsigned int J=0; J<DIM; ++J){
      T e1_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J)+F[i][2]*(2==J))/sqrt(3.0);
      T e2_FiJ=(F[i][0]*(0==J)-F[i][1]*(1==J))/sqrt(2.0);
      T e3_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))/sqrt(6.0);
      T e4_FiJ=(F[i][2]*(1==J)+F[i][1]*(2==J))/2.0;
      T e5_FiJ=(F[i][0]*(2==J)+F[i][2]*(0==J))/2.0;
      T e6_FiJ=(F[i][1]*(0==J)+F[i][0]*(1==J))/2.0;
      T e2_1_FiJ=((0==J)*dF[i][0][0]-(1==J)*dF[i][1][0])/sqrt(2.0);
      T e2_2_FiJ=((0==J)*dF[i][0][1]-(1==J)*dF[i][1][1])/sqrt(2.0);
      T e3_1_FiJ=((0==J)*dF[i][0][0]+(1==J)*dF[i][1][0]-2*(2==J)*dF[i][2][0])/sqrt(6.0);
      T e3_2_FiJ=((0==J)*dF[i][0][1]+(1==J)*dF[i][1][1]-2*(2==J)*dF[i][2][1])/sqrt(6.0);

      //P
      P[i][J]=(2*Eii*e1)*e1_FiJ						\
	+(2*E2c*e2 -6*E3c*e2*e3 + 4*E4c*e2*(e2*e2+e3*e3))*e2_FiJ		\
	+(2*E2c*e3 +3*E3c*(e3*e3-e2*e2) + 4*E4c*e3*(e2*e2+e3*e3))*e3_FiJ	\
	+(2*Eij*e4)*e4_FiJ						\
	+(2*Eij*e5)*e5_FiJ						\
	+(2*Eij*e6)*e6_FiJ						\
	+(2*Eg*e2_1+(2.0/sqrt(3.0))*Eg*e3_1 - sqrt(3.0)*Eh*e3_1)*e2_1_FiJ \
	+(2*Eg*e2_2-(2.0/sqrt(3.0))*Eg*e3_2 + sqrt(3.0)*Eh*e3_2)*e2_2_FiJ \
	+(2*Eg*e3_1/3.0 + (2.0/sqrt(3.0))*Eg*e2_1 + 2*Eh*e3_1 - sqrt(3)*Eh*e2_1)*e3_1_FiJ \
	+(2*Eg*e3_2/3.0 - (2.0/sqrt(3.0))*Eg*e2_2 + 2*Eh*e3_2 + sqrt(3)*Eh*e2_2)*e3_2_FiJ;
	
      //gradient terms
      for (unsigned int K=0; K<DIM; ++K){
	T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K)/sqrt(2.0);
	T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K)/sqrt(2.0);
	T e3_1_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(0==K)/sqrt(6.0);
	T e3_2_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(1==K)/sqrt(6.0);
	//Beta
	Beta[i][J][K]= (2*Eg*e2_1 + 2*Eg*e3_1/sqrt(3.0) - sqrt(3.0)*Eh*e3_1)*e2_1_FiJK \
	  + (2*Eg*e2_2 - 2*Eg*e3_2/sqrt(3.0) + sqrt(3.0)*Eh*e3_2)*e2_2_FiJK \
	  + (2*Eg*e3_1/3.0 + 2*Eg*e2_1/sqrt(3.0) + 2*Eh*e3_1 - sqrt(3.0)*Eh*e2_1)*e3_1_FiJK \
	  + (2*Eg*e3_2/3.0 - 2*Eg*e2_2/sqrt(3.0) + 2*Eh*e3_2 + sqrt(3.0)*Eh*e2_2)*e3_2_FiJK;
      }
    }
  }
  
  //chemical potential
  T mu = 4*C4*c*c*c+2*C2*c			\
    + (E4/(2*Cs))*(e2*e2+e3*e3)*(e2*e2+e3*e3)	\
    + (E3/(2*Cs))*e3*(e3*e3-3*e2*e2)		\
    + (-E2/(2*Cs))*(e2*e2+e3*e3);
  T dmuc= 12*C4*c*c + 2*C2;
  T dmue2= (E4/(2*Cs))*4*(e2*e2+e3*e3)*e2 + (E3/(2*Cs))*e3*(-6*e2) + (-E2/(2*Cs))*(2*e2);
  T dmue3= (E4/(2*Cs))*4*(e2*e2+e3*e3)*e3 + (E3/(2*Cs))*(3*e3*e3 - 3*e2*e2) + (-E2/(2*Cs))*(2*e3);
  
  /* //get shape function values */
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];

  //Compute Residual
  bool surfaceFlag=p->atboundary;
  T (*Ra)[DIM+1] = (T (*)[DIM+1])R;
  for (unsigned int a=0; a<(unsigned int)nen; a++) {
    double N1[DIM], N2[DIM][DIM];
    for (unsigned int i=0; i<DIM; i++){
      N1[i]=Nx[a][i];
      for (unsigned int j=0; j<DIM; j++){
	N2[i][j]=Nxx[a][i][j];
      }
    }

    //Chemistry
    T laplace_c=0;
    for (unsigned int i=0; i<DIM; i++) laplace_c+=cxx[i][i];
    T Rc=0.0;
    if (!surfaceFlag){
      // Na * c_t
      Rc += N[a] * (c-c0)*(1.0/user->dt);
      // grad(Na) . D*(dmuc*grad(C)+dmue2*grad(e2)+dmue3*grad(e3))
      double laplace_N=0.0;
      for (unsigned int i=0; i<DIM; i++){
	T e2x=0, e3x=0;
	switch (i) {
	case 0:
	  e2x=e2_1; e3x=e3_1; break;
	case 1:
	  e2x=e2_2; e3x=e3_2; break;
	case 2:
	  e2x=e2_3; e3x=e3_3; break;
	} 
	Rc += N1[i]*D*(dmuc*cx[i]+dmue2*e2x+dmue3*e3x);
	laplace_N += N2[i][i];
      }
      // lambda * del2(Na) * D * del2(c)
      Rc += CLambda*laplace_N*D*laplace_c;
    }
    else{
      // -grad(Na) . (D*del2(c)) n
      T t1 = D*laplace_c;
      double laplace_N=0.0;
      for (unsigned int i=0; i<DIM; i++){
	Rc += -N1[i]*t1*n[i];
	laplace_N += N2[i][i];
      }
      // -(gamma*del2(Na)*D)*grad(C).n
      T t2 = gamma*laplace_N*D;
      for (unsigned int i=0; i<DIM; i++){
	Rc += -t2*cx[i]*n[i];
      }
      // (C/he)*(grad(Na).n)*D*(grad(C).n)
      double t3=0.0;
      T t4 = (C/he)*D;
      T t5=0.0;
      for (unsigned int i=0; i<DIM; i++){
	t3 += N1[i]*n[i];
	t5 += cx[i]*n[i];
      }
      Rc += t3*t4*t5;
      Rc *= CLambda;
      //flux term
      // Na*J
      Rc += -N[a]*flux;
    }
    Ra[a][0] = Rc;

    //Mechanics
    if (!surfaceFlag) {
      for (unsigned int i=0; i<DIM; i++){
	T Ru_i=0.0;
	for (unsigned int j=0; j<DIM; j++){
	  //grad(Na)*P
	  Ru_i += N1[j]*P[i][j];
	  for (unsigned int k=0; k<DIM; k++){
	    Ru_i += N2[j][k]*Beta[i][j][k];
	  }
	}
	Ra[a][i+1] = Ru_i;
      }
    }
  }
  return 0;
}


#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint p,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx)
{
  Function(p, dt, shift, V, t, U, t0, U0, R, ctx);
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "Jacobian"
PetscErrorCode Jacobian(IGAPoint p,PetscReal dt,
				PetscReal shift,const PetscScalar *V,
				PetscReal t,const PetscScalar *U,
				PetscReal t0,const PetscScalar *U0,
				PetscScalar *K,void *ctx)
{
  const PetscInt nen=p->nen, dof=DIM+1;
  doubleAD *U_AD=new doubleAD[nen*(DIM+1)];
  const PetscReal (*U2)[DIM+1] = (PetscReal (*)[DIM+1])U;
  for(int n=0; n<nen; n++){
    for(int d=0; d<dof; d++){
      U_AD[n*dof+d]=U2[n][d];
      U_AD[n*dof+d].diff(n*dof+d, dof*nen);
    }
  }
 
  doubleAD *R= new doubleAD[nen*dof];
  Function (p, dt, shift, V, t, U_AD, t0, U0, R, ctx);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = R[n1*dof+d1].dx(n2*dof+d2);
	}
      }
    }				
  }
  delete []R;
  delete []U_AD;
  return 0;    
}


PetscErrorCode L2Residual(IGAPoint p,const PetscScalar *U, PetscScalar *R, void *ctx)
{	
  AppCtx *user = (AppCtx *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  //displacement field variables
  double u[DIM], ux[DIM][DIM];
  computeField<double,DIM,DIM+1>(VECTOR,1,p,U,&u[0],&ux[0][0]);

  //Compute F
  PetscReal F[3][3];
  for (PetscInt i=0; i<3; i++) {
    for (PetscInt j=0; j<3; j++) {
      if ((i<DIM) && (j<DIM))  F[i][j]=(i==j)+ux[i][j];
      else  F[i][j]=(i==j);
    }
  }
  
  //Compute strain metric, E  (E=0.5*(F^T*F-I))
  PetscReal E[3][3];
  for (PetscInt i=0; i<3; i++){
    for (PetscInt j=0; j<3; j++){
      E[i][j] = -0.5*(i==j);
      for (PetscInt k=0; k<3; k++){
	E[i][j] += 0.5*F[k][i]*F[k][j];
      }
    }
  }
  //new strain metrics
  PetscReal e1=(E[0][0]+E[1][1]+E[2][2])/sqrt(3.0);
  PetscReal e2=(E[0][0]-E[1][1])/sqrt(2.0);
  PetscReal e3=(E[0][0]+E[1][1]-2*E[2][2])/sqrt(6.0);
  PetscReal e4=E[1][2], e5=E[2][0], e6=E[0][1];

  //compute distance to nearest well
  PetscReal Es=user->Es;
  PetscReal x[3],y[3]; 
  x[0]=0; y[0]=Es; //first well 
  x[1]=-Es*cos(30.0*PI/180.0); y[1]=-Es*sin(30.0*PI/180.0); //second well
  x[2]=+Es*cos(30.0*PI/180.0); y[2]=-Es*sin(30.0*PI/180.0); //third well
  PetscReal dist=sqrt(std::pow(e2-x[0],2.0)+std::pow(e3-y[0],2.0));
  unsigned int wellID=1; 
  for(unsigned int i=1; i<3; i++){
    if(dist>sqrt(pow(e2-x[i],2.0)+pow(e3-y[i],2.0))){
      dist=sqrt(pow(e2-x[i],2.0)+pow(e3-y[i],2.0));
      wellID=i+1;
    }
  }

  //store L2 projection residual
  const PetscReal *N;
  IGAPointGetShapeFuns(p,0,(const PetscReal**)&N);
  PetscReal (*RLocal)[DIM+1] = (PetscReal (*)[DIM+1])R;
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      PetscReal val=0.0;
      switch (d1) {
      case 0:
	val=e2; break;
      case 1:
	val=e3; break;
      case 2:
	val=dist; break;
      case 3: //only in 3D
	val=(PetscReal) wellID; break;
      } 
      RLocal[n1][d1] = N[n1]*val;
    }
  }
  return 0;
}

PetscErrorCode L2Jacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx)
{	
  AppCtx *user = (AppCtx *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //store L2 projection matrix
  const PetscReal *N;
  IGAPointGetShapeFuns(p,0,(const PetscReal**)&N);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
	  PetscScalar val=0.0;
	  if (d1==d2) {val = N[n1] * N[n2];}
	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = val;
	}
      }
    }
  }
  return 0;
}


PetscErrorCode WriteSolution(Vec U, int step, IGA iga, const char filePrefix[], bool restartFile=false)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  char	filename[256];
  sprintf(filename,"%s-%d.dat",filePrefix,step);
  ierr = IGAWriteVec(iga, U, filename);CHKERRQ(ierr);
  //restart file
  if (restartFile){
     MPI_Comm        comm;
     PetscViewer     viewer;
     sprintf(filename,"res%s-%d.dat",filePrefix,step);
     ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
     ierr = PetscViewerBinaryOpen(comm,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
     ierr = VecView(U,viewer);CHKERRQ(ierr);
     ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
     ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ProjectSolution(IGA iga, Vec U, PetscInt step, const char filePrefix[], void *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //Setup linear system for L2 Projection
  //Ideally we would use KSP here as this is a linear problem... \\
  //but since we need to pass solution vector U we are forced to use SNES context instead.
  
  // Setup the nonlinear solver
  SNES snes;
  Vec Utemp;;
  ierr = IGASetFormFunction(iga,L2Residual,user);CHKERRQ(ierr);
  ierr = IGASetFormJacobian(iga,L2Jacobian,user);CHKERRQ(ierr);
  ierr = IGACreateSNES(iga,&snes);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&Utemp);CHKERRQ(ierr);
  ierr = VecZeroEntries(Utemp);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,Utemp);CHKERRQ(ierr);
  //write solution
  char filename[256];
  sprintf(filename,"out%sE",filePrefix);
  ierr = WriteSolution(Utemp, step, iga, filename);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = VecDestroy(&Utemp);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
}

typedef struct {
  PetscReal c, ux, uy;
#if DIM==3  
  PetscReal uz;
#endif
} Field;

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,1+DIM,&da);CHKERRQ(ierr);
#if DIM==3
  Field ***u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscInt i,j,k;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
	u[k][j][i].c = user->cbar + 0.01*(0.5 - (double)(std::rand() % 100 )/100.0);
	u[k][j][i].ux=0.0;
	u[k][j][i].uy=0.0;
	u[k][j][i].uz=0.0;
      }
    }
  }
#elif DIM==2
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      u[j][i].c = user->cbar + 0.01*(0.5 - (double)(std::rand() % 100 )/100.0);
      u[j][i].ux=0.0;
      u[j][i].uy=0.0;
    }
  }
#endif
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "OutputMonitor"
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  sprintf(filename,"./out%d.dat",it_number);
  ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
  double startTime = MPI_Wtime();
 
  /* Define simulation specific parameters */
  AppCtx user;
  
  //problem parameters
  user.Cs=1.0;
  user.Cd=1.0;
  user.C4=user.Cd/pow(user.Cs,4.0); 
  user.C2=-2*user.Cd/pow(user.Cs,2.0);
  user.CLambda=1.0e-2;
  user.Es=0.1;
  user.Ed=1.0;
  user.E4=3*user.Ed/pow(user.Es,4.0);
  user.E3=2*user.Ed/pow(user.Es,3.0)-2*user.E4*user.Es;
  user.E2=user.Ed/pow(user.Es,2.0);
  user.Eii=user.Ed/pow(user.Es,2.0);
  user.Eij=user.Ed/pow(user.Es,2.0);
  user.El=0.0;
  user.Eg=pow(user.El,2.0)*user.Ed/pow(user.Es,2.0);
  user.Eh=pow(user.El,2.0)*user.Ed/pow(user.Es,2.0);
  user.D=1.0;
  user.flux=1.0;
  user.cbar=0.0;
  user.gamma=1.0;
  user.C=5.0;
  user.dt=1.0e-5;

  /* Set discretization options */
  PetscInt nsteps = 100;
  PetscInt N=10, p=2, C=PETSC_DECIDE, resStep=0;
  PetscBool output = PETSC_FALSE; 
  PetscBool monitor = PETSC_FALSE; 
  char filePrefix[PETSC_MAX_PATH_LEN] = {0};
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","CahnHilliard2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-N","number of elements (along one DIMension)",__FILE__,N,&N,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-p","polynomial order",__FILE__,p,&p,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-C","global continuity order",__FILE__,C,&C,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-file_prefix","File Prefix",__FILE__,filePrefix,filePrefix,sizeof(filePrefix),PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-ch_output","Enable output files",__FILE__,output,&output,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-ch_monitor","Compute and show statistics of solution",__FILE__,monitor,&monitor,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nsteps","Number of load steps to take",__FILE__,nsteps,&nsteps,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-dt","time step",__FILE__,user.dt,&user.dt,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-res_step","Restart Step",__FILE__,resStep,&resStep,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;
  //
  user.he=1.0/N;
  //PetscPrintf(PETSC_COMM_WORLD,"\nLambda_u value is: %8.2e\n",user.lambda_u);

  //
  if (p < 2 || C < 0) /* Problem requires a p>=2 C1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Problem requires minimum of p = 2");
  if (p <= C)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(iga,1+DIM);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,C);CHKERRQ(ierr);
  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisCopy(axis0,axis1);CHKERRQ(ierr);
#if DIM==3
  IGAAxis axis2;
  ierr = IGAGetAxis(iga,2,&axis2);CHKERRQ(ierr);
  ierr = IGAAxisCopy(axis0,axis2);CHKERRQ(ierr);
#endif
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;
  if (resStep==0){
    char meshfilename[256];
    sprintf(meshfilename, "mesh.dat");
    PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
    ierr = IGAWrite(iga, meshfilename);CHKERRQ(ierr);
  }
  
  Vec U,U0;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U0);CHKERRQ(ierr);
  if (resStep>0){
    MPI_Comm comm;
    PetscViewer viewer;
    char restartfilename[256];
    ierr = PetscObjectGetComm((PetscObject)U0,&comm);CHKERRQ(ierr);
    sprintf(restartfilename,"res%s-%d.dat",filePrefix,resStep);
    ierr = PetscViewerBinaryOpen(comm,restartfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(U0,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    PetscPrintf(PETSC_COMM_WORLD,"\nReading solution from restart file: %s\n",restartfilename);
  }
  else{
    ierr = FormInitialCondition(iga, U0, &user); 
  }
  ierr = VecCopy(U0, U);CHKERRQ(ierr);

  //Dirichlet BC
  double dVal=0.01*user.Es;
  ierr = IGASetBoundaryValue(iga,0,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,0,2,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,1,dVal);CHKERRQ(ierr);
#if DIM==3
  ierr = IGASetBoundaryValue(iga,2,0,3,0.0);CHKERRQ(ierr);
#endif
  //
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);
  //
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,10000,10.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  if (output) {
    ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  }
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

#if PETSC_VERSION_LE(3,3,0)
  ierr = TSSolve(ts,U,NULL);CHKERRQ(ierr);
#else
  ierr = TSSolve(ts,U);CHKERRQ(ierr);
#endif
  //
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

