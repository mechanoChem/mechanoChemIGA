#include "genericHeaders.h"
//extern "C" {
#include "petiga.h"
//}
#include "../src/configurationalForces/physicsHeaders.h"
#include "../applications/configurationalForces/3D/parameters.h"

//include automatic differentiation library
#ifdef ADSacado
#include <Sacado.hpp>
#else
#include "adept.h"
#endif

#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint p,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx)
{
  Function<PetscReal>(p, dt, shift, V, t, U, t0, U0, R, ctx);
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "Jacobian"
template<unsigned int DOF>
PetscErrorCode Jacobian(IGAPoint p,PetscReal dt,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscScalar *U,
			PetscReal t0,const PetscScalar *U0,
			PetscScalar *K,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;
  const PetscInt nen=p->nen, dof=DOF;
  const PetscReal (*U2)[DOF] = (PetscReal (*)[DOF])U;
#ifdef ADSacado

  typedef Sacado::Fad::SFad<double,numVars> doubleAD;
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
  typedef adept::adouble doubleAD;
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

template PetscErrorCode Jacobian<3>(IGAPoint p,PetscReal dt,
				    PetscReal shift,const PetscScalar *V,
				    PetscReal t,const PetscScalar *U,
				    PetscReal t0,const PetscScalar *U0,
				    PetscScalar *K,void *ctx);
template PetscErrorCode Jacobian<4>(IGAPoint p,PetscReal dt,
				    PetscReal shift,const PetscScalar *V,
				    PetscReal t,const PetscScalar *U,
				    PetscReal t0,const PetscScalar *U0,
				    PetscScalar *K,void *ctx);
template PetscErrorCode Jacobian<6>(IGAPoint p,PetscReal dt,
				    PetscReal shift,const PetscScalar *V,
				    PetscReal t,const PetscScalar *U,
				    PetscReal t0,const PetscScalar *U0,
				    PetscScalar *K,void *ctx);
