/**
 * @page example6 Example 6 : Unconditionally stable scheme for gradient elasticity
 * \dontinclude gradElastTime/userFunctions.cc
 * As in previous examples, we include the header file declaring the required user functions.
 *
 * \line userFunction
 *
 * We define the initial conditions.
 *
 * <b> The \c uinit function </b>
 *
 * The initial conditions for the vector field are defined for this example problem as following.  
 *
 * \skip template
 * \until //end
 *
 * <b> The \c defineParameters function </b>
 *
 * We are going to use periodic boundary conditions, which can be set as:
 *
 * \skip user.periodic[0]
 * \until user.periodic[2]
 *
 *
 * Here, we define the mesh by setting the number of elements in each direction, e.g. a 4x4x4 element mesh.
 *
 * \skip user.N[0]
 * \until user.N[2]
 *
 * We also define the dimensions of the domain, e.g. a unit cube.
 *
 * \skip user.L[0]
 * \until user.L[2]
 *
 * We specify the number of vector and scalar solution by adding the name of each field to
 * their respective vector. Here, we only have one vector solution field (the displacement).
 *
 * \skip "displacement"
 * \until "displacement"
 * 
 * We also specify the polynomial order of the basis splines and the global continuity.
 * 
 * \skip polyOrder
 * \until globalContinuity
 *
 * We redirect the desired user function pointer to the \c uinit functions that we
 * defined above.
 *
 * \skip InitialConditions
 * \until = uinit
 *
 * Finally, we define various (9) material parameters that describe our gradient elasticity.
 *
 * \skip malloc
 * \until par_mat[8]
 *
 * <b> The \c residual function </b>
 *
 * The residual function for an unconditionally stable second-order scheme for gradient-elasticity is used in this example.
 *
 * We first declare \c "eval_residual" non-member function to be used in the member function, \c "residual".
 * 
 * \skip template
 * \until eval_residual
 *
 * The definition of the \c eval_residual function is postponed until the end of the file as it is lengthy.
 * 
 * For a complex problem like this example, where  it is convenient to unfold "u.val", "u.grad", and "u.hess" and put them into a single array, ui[]. 
 *
 * \skip ui[]=
 * \until };
 * 
 * We do the same for previous solutions represented by ".valP", ".gradP", and ".hessP", and ".valPP", ".gradPP", and ".hessPP" as well as for the test functions "w1" and "w2" and produce arrays, up[], upp[], and w[], respectively.
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
 * A snippet of the code
 * ==============================
 * 
 * \dontinclude gradElastTime/userFunctions.cc
 * \skip userFunctions.h
 * \until comp+63
 *
 *
 * 
 */
