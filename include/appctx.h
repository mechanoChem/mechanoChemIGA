#ifndef appctx_
#define appctx_

typedef struct {
  IGA iga;
  Vec *U, *U0;
  PetscReal dt;
  PetscReal he;
  AppCtxKSP* appCtxKSP;
  PetscReal f0Norm;
} AppCtx;

#endif
