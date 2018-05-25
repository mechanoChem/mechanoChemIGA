/**
 * @page example6 Example 6 : Gradient elasticity with periodic boundary conditions.
 * \dontinclude gradElast/userFunctions.cc
 * As in previous examples, we include the header file declaring the required user functions.
 *
 * \line userFunction
 *
 * Then, we define the initial condition.
 *
 * <b> The \c uinit function </b>
 *
 * The initial condition for the vector field is defined for this example problem to introduce small perturbation as following.  
 *
 * \skip template
 * \until //end
 *
 * <b> The \c defineParameters function </b>
 * 
 * Instead of applying Dirichlet boundary conditions, in this example, we apply periodic conditions as:
 *
 * \skip user.periodic[0]
 * \until user.periodic[2]
 * 
 * We also define the mesh as in other examples by setting the number of elements in each direction.
 *
 * \skip user.N[0]
 * \until user.N[2]
 *
 * We also define the dimensions of the domain, e.g. a unit cube.
 *
 * \skip user.L[0]
 * \until user.L[2]
 *
 * We specify the number of vector solution by adding the name of the field to a vector. 
 *
 * \line "displacement"
 * 
 * We also specify the polynomial order of the basis splines and the global continuity.
 * 
 * \skip polyOrder
 * \until globalContinuity
 *
 * We redirect the desired user function pointer \c uinit function that we defined above.
 *
 * \line InitialConditions
 *
 * Finally, we define various (9) material parameters that describe the gradient elasticity.
 *
 * \skip malloc
 * \until par_mat[8]
 *
 * <b> The \c residual function </b>
 *
 * The residual function for the gradient elasticity with periodic conditions is used in this example.
 *
 * We first declare \c "eval_residual" non-member function to be used in the member function, \c "residual".
 * 
 * \skip template
 * \until eval_residual
 *
 * The definition of the \c eval_residual function is postponed until the end of the file as it is lengthy.
 * 
 * It is convenient to unfold "u.val", "u.grad", and "u.hess" and put them into a single array, ui[]. 
 *
 * \skip ui[]=
 * \until };
 * 
 * We do the same for previous solutions represented by ".valP", ".gradP", and ".hessP" as well as for the test functions "w2" and produce arrays, u0[] and w[], respectively.
 * 
 * We then evaluate the residual vector at a given quadrature point (residual[]) using the declared function "eval_residual".
 *
 * \skip T residual[
 * \until eval
 *
 * Finally, we multiply residual[] by test functions and form the integrand of the weak form at the given quadrature point. 
 *
 * \skip r=0
 * \until dV
 *
 * The complete code
 * ==============================
 * 
 * \dontinclude gradElast/userFunctions.cc
 * \skip userFunctions.h
 * \until "userFunctionsInstantiation.h"
 *
 *
 * 
 */
