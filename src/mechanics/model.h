#ifndef mechanicsND_
#define mechanicsND_

//residual function implementation
#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;
  PetscReal c=user->lambda; //used as a load paramter 0<c<1

  //displacement field variable
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DIM>(VECTOR,0,p,U,&u[0],&ux[0][0],&uxx[0][0][0]);

  //Compute F (I+Ux), dF (Uxx)
  T F[DIM][DIM], dF[DIM][DIM][DIM];
  for (unsigned int i=0; i<DIM; i++) {
    for (unsigned int J=0; J<DIM; J++) {
      F[i][J]=(i==J)+ux[i][J];
      for (unsigned int K=0; K<DIM; K++) {
	dF[i][J][K]=uxx[i][J][K];
      }
    }
  }

  //Compute strain metric, E=0.5*(F^T*F-I)
  T E[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
    }
  }

  //compute P and Beta
  T P[DIM][DIM], Beta[DIM][DIM][DIM];
 
  //2D model
#if DIM==2  
  //new strain metrics (2D)
  T e1=(E[0][0]+E[1][1]);
  T e2=(E[0][0]-E[1][1]);
  T e6=(E[0][1]);
  T e2_1=0.0, e2_2=0.0; 
  for (unsigned int i=0; i<DIM; ++i){
    e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0]);
    e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1]);
  }
  //
  for (unsigned int i=0; i<DIM; ++i){
    for (unsigned int J=0; J<DIM; ++J){
      T e1_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J));
      T e2_FiJ=(F[i][0]*(0==J)-F[i][1]*(1==J));
      T e6_FiJ=(F[i][1]*(0==J)+F[i][0]*(1==J))/2.0;
      T e2_1_FiJ=((0==J)*dF[i][0][0]-(1==J)*dF[i][1][0]);
      T e2_2_FiJ=((0==J)*dF[i][0][1]-(1==J)*dF[i][1][1]);
      //P
      P[i][J]=PiJ;

      //gradient terms
      for (unsigned int K=0; K<DIM; ++K){
	T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K);
	T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K);
	//Beta
	Beta[i][J][K]=BetaiJK;
      }
    }
  }
  
  //3D model
#elif DIM==3
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
  //
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
      T e2_3_FiJ=((0==J)*dF[i][0][2]-(1==J)*dF[i][1][2])/sqrt(2.0);     
      T e3_1_FiJ=((0==J)*dF[i][0][0]+(1==J)*dF[i][1][0]-2*(2==J)*dF[i][2][0])/sqrt(6.0);
      T e3_2_FiJ=((0==J)*dF[i][0][1]+(1==J)*dF[i][1][1]-2*(2==J)*dF[i][2][1])/sqrt(6.0);
      T e3_3_FiJ=((0==J)*dF[i][0][2]+(1==J)*dF[i][1][2]-2*(2==J)*dF[i][2][2])/sqrt(6.0);
      //P
      P[i][J]=PiJ;   
      
      //gradient terms
      for (unsigned int K=0; K<DIM; ++K){
	T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K)/sqrt(2.0);
	T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K)/sqrt(2.0);
	T e2_3_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(2==K)/sqrt(2.0);
	T e3_1_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(0==K)/sqrt(6.0);
	T e3_2_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(1==K)/sqrt(6.0);
	T e3_3_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(2==K)/sqrt(6.0);
	//Beta
	Beta[i][J][K]=BetaiJK;
      }
    }
  }
#else
  PetscPrintf(PETSC_COMM_WORLD,"only material models for DIM=2, DIM=3 implemented.... but DIM input is %u\n",DIM); 
      exit(-1);      	
#endif   

  //get shape function values
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];

  //compute residuals
  bool surfaceFlag=p->atboundary;
  for (unsigned int a=0; a<(unsigned int)nen; a++) {
    double N1[DIM], N2[DIM][DIM];
    for (unsigned int i=0; i<DIM; i++){
      N1[i]=Nx[a][i];
      for (unsigned int j=0; j<DIM; j++){
	N2[i][j]=Nxx[a][i][j];
      }
    }

    //mechanics
    if (!surfaceFlag) {
      for (unsigned int i=0; i<DIM; i++){
	T Ru_i=0.0;
	for (unsigned int J=0; J<DIM; J++){
	  //grad(Na)*P
	  Ru_i += N1[J]*P[i][J];
	  for (unsigned int K=0; K<DIM; K++){
	    Ru_i += N2[J][K]*Beta[i][J][K];
	  }
	}
	R[a*dof+i] = Ru_i;
      }
    }
  }
  return 0;
}

#endif
