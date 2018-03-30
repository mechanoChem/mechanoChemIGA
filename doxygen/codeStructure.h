/**
 * @page codestructure Code Structure
 *
 * The code is divided into two parts: user defined functions and core functions.
 * The user defined functions is a list of functions that can/must be defined by the user.
 * The \c defineParameters and \c residual functions (italicized) require a user definition;
 * the remaining functions have default definitions that can be overridden by the user if desired.
 * These user functions are called by the core functions, as shown in the following diagram.
 *
 * @dot
 digraph G
{
  graph[rankdir="LR",bgcolor="transparent",splines="ortho"];

  subgraph cluster_userFunctions{
   node [fontname="FreeSans",fontsize=13,fontcolor="blue",
        shape=record,height=0.2,width=0.4,
        color="black", fillcolor="white", style="filled"];
   fontname="FreeSans";
   fontsize=15;
   label=< <u>  User level 0  </u><br/><br/>User defined functions>;
   color="transparent";
   scalarIC [label="scalarInitialConditions",URL="\ref scalarInitialConditions"];
   vectorIC [label="vectorInitialConditions",URL="\ref vectorInitialConditions"];
   defineParameters [label=<{<B><I>defineParameters</I></B>}>,URL="\ref defineParameters"];
   BCs [label="boundaryConditions",URL="\ref boundaryConditions"];
   residual [label=<{<B><I>residual</I></B>}>,URL="\ref residual"];
   loadStep [label="loadStep",URL="\ref loadStep"];
   adaptiveTS [label="adaptiveTimeStep",URL="\ref adaptiveTimeStep"];
   projectFields [label="projectFields",URL="\ref projectFields"];
   {rank=same; scalarIC, vectorIC, BCs, loadStep, adaptiveTS, residual, defineParameters, projectFields}
  }

  subgraph cluster_baseFunctions{
   node [fontname="FreeSans",fontsize=13,fontcolor="blue",
        shape=record,height=0.2,width=0.4,
        color="black", fillcolor="white", style="filled"];
   fontname="FreeSans";
   fontsize=15;
   label = < <u>  Developer level  </u><br/><br/>Core functions>;
   color="transparent";
   Space1 [color="transparent",label=""];
   formIC [label="FormInitialCondition",URL="\ref FormInitialCondition"];
   setup [label="Setup",URL="\ref Setup"];
   Space2 [color="transparent",label=""];
   quadPtResidual [label="QuadPtResidual",URL="\ref QuadPtResidual"];
   StepUpdate [label="StepUpdate",URL="\ref StepUpdate"];
   Space3 [color="transparent",label=""];
   ProjectionResidual [label="ProjectionResidual",URL="\ref ProjectionResidual"];
   {rank=same; Space1, formIC, setup, Space2, quadPtResidual, StepUpdate, Space3, ProjectionResidual}
  }

  subgraph cluster_worflow{
   node [fontname="FreeSans",fontsize=13,
        color="black", fillcolor="white", style="filled"];
   fontname="FreeSans";
   fontsize=15;
   label = <Workflow<br/><font color="white">.</font> >;
   color="transparent";
   SETUP [label="\nSETUP\n\n"];
   ASSMBLY [label="\nASSEMBLY\n\n"];
   SOLVE [label="\nSOLVE\n\n"];
   PP [label="\nPOSTPROCESSING\n\n"];
   {rank=same; SETUP, ASSMBLY, SOLVE, PP}
  }

  d1 [shape=point,width=0.001,height=0.001];
  d2 [shape=point,width=0.001,height=0.001];
  d3 [shape=point,width=0.001,height=0.001];

  {scalarIC,vectorIC} -> d1 [color="black",fontsize=10,style="solid",fontname="FreeSans",dir=none];
  d1 -> formIC [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  formIC -> setup [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  {defineParameters,BCs} -> d2 [color="black",fontsize=10,style="solid",fontname="FreeSans",dir=none];
  d2 -> setup [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  residual -> quadPtResidual [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  {loadStep,adaptiveTS} -> d3 [color="black",fontsize=10,style="solid",fontname="FreeSans",dir=none];
  d3 -> StepUpdate [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  projectFields -> ProjectionResidual [color="black",fontsize=10,style="solid",fontname="FreeSans"];
  formIC -> SETUP [color="transparent",fontsize=10,style="dashed",arrowhead="none",fontname="FreeSans"];
  quadPtResidual -> ASSMBLY [color="transparent",fontsize=10,style="solid",fontname="FreeSans"];
  StepUpdate -> SOLVE [color="transparent",fontsize=10,style="solid",fontname="FreeSans"];
  ProjectionResidual -> PP [color="transparent",fontsize=10,style="solid",fontname="FreeSans"];
  SETUP -> ASSMBLY [color="black",fontsize=10,style="solid",fontname="FreeSans",dir=none];
  ASSMBLY -> SOLVE [color="black",fontsize=10,style="solid",fontname="FreeSans",dir=none];
  SOLVE -> PP [color="black",fontsize=10,style="solid",fontname="FreeSans",dir=none];
}
 * @enddot
 *
 */
