#ifndef appctx_
#define appctx_

typedef struct {
  IGA iga;
  TS* ts;
  Vec *U, *U0;
  PetscReal dt;
  PetscReal he;
  PetscReal norm;
} AppCtx;

#endif
