//include automatic differentiation library
#ifdef ADSacado
#include <Sacado.hpp>
typedef Sacado::Fad::SFad<double,numVars> doubleAD;
//typedef Sacado::Fad::DFad<double> doubleAD;
#else
#include "adept.h"
typedef adept::adouble doubleAD;
#endif

template <class T>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx);

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
  AppCtx *user = (AppCtx *)ctx;
  const PetscInt nen=p->nen, dof=DIM+1;
  const PetscReal (*U2)[DIM+1] = (PetscReal (*)[DIM+1])U;
#ifdef ADSacado
  if (dof*nen!=numVars) {
    PetscPrintf(PETSC_COMM_WORLD,"\ndof*nen!=numVars.... Set numVars = %u\n",dof*nen); exit(-1);
  }
  std::vector<doubleAD> U_AD(nen*dof);
  for(int i=0; i<nen*dof; i++){
    U_AD[i]=U[i];
    U_AD[i].diff(i, dof*nen);
  } 
  std::vector<doubleAD> R(nen*dof);
  Function<doubleAD> (p, dt, shift, V, t, &U_AD[0], t0, U0, &R[0], ctx);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
      	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = R[n1*dof+d1].dx(n2*dof+d2);
	}
      }
    }				
  }
#else
  adept::Stack s;
  std::vector<doubleAD> U_AD(nen*dof);
  adept::set_values(&U_AD[0],nen*dof,U);
  s.new_recording();
  std::vector<doubleAD> R(nen*dof);
  Function<doubleAD> (p, dt, shift, V, t, &U_AD[0], t0, U0, &R[0], ctx);
  s.independent(&U_AD[0],nen*dof);
  s.dependent(&R[0],nen*dof);
  std::vector<double> K2(nen*dof*nen*dof);
  s.jacobian(&K2[0]);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
      	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] =  K2[n2*dof*nen*dof + d2*nen*dof + n1*dof + d1];
	}
      }
    }				
  }
#endif  
  return 0;    
}

