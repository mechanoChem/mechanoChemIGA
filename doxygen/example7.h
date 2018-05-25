/**
 * @page example7 Example 7 : Unconditionally stable scheme for mechano-chemistry
 * \dontinclude mechanoChemStab/userFunctions.cc
 * As in previous examples, we include the header file declaring the required user functions.
 *
 * \line userFunction
 *
 * Then, we define initial/boundary conditions.
 *
 * <b> The \c boundaryConditions function </b>
 *
 * We apply boundary conditions so that the material is free to move in tangential directions except on z=1, where we have the free-surface condition.  
 *
 * \skip boundary
 * \until //end
 *
 * <b> The \c cmuinit function and the \c uinit function </b>
 *
 * The initial conditions for the scalar and the vector fields are defined for this example problem as following.  
 *
 * \skip template
 * \until //end
 *
 * \skip template
 * \until //end
 *
 * <b> The \c defineParameters function </b>
 *
 * Here, we define the mesh by setting the number of elements in each direction, e.g. a 64x64x64 element mesh.
 *
 * \skip user.N[0]
 * \until user.N[2]
 *
 * We also define the dimensions of the domain, e.g. a unit cube.
 *
 * \skip user.L[0]
 * \until user.L[2]
 *
 * We specify the number of vector and scalar solution and projection fields by adding the name of each field to
 * their respective vector. Here, we have one vector solution field (the displacement) and two scalar projection field
 * (the chemical composition and the chemical potential).
 *
 * \skip "c"
 * \until "displacement"
 * 
 * We also specify the polynomial order of the basis splines and the global continuity.
 * 
 * \skip polyOrder
 * \until globalContinuity
 *
 * We redirect the desired user function pointers to the \c boundaryConditions, \c cmuinit, and \c uinit functions that we
 * defined above.
 *
 * \skip boundaryConditions
 * \until = uinit
 *
 * We then define the timestep size.
 *
 * \skip 1.e-4
 * \line dtVal
 *
 * Finally, we define various (25) material parameters that describe the mechano-chemistry.
 *
 * \skip malloc
 * \until par_mat[24]
 *
 * <b> The \c residual function </b>
 *
 * The residual function for an unconditionally stable second-order scheme for mechano-chemistry is used in this example.
 *
 * We first declare \c "eval_residual" non-member function to be used in the member function, \c "residual".
 * 
 * \skip template
 * \until eval_residual
 *
 * The definition of the \c eval_residual function is postponed until the end of the file as it is lengthy.
 * 
 * For a complex problem like this example, where  it is convenient to unfold "c.val" and "u.val", "c.grad" and "u.grad", and "c.hess" and "u.hess" and put them into a single array, ui[]. 
 *
 * \skip ui[]=
 * \until };
 * 
 * We do the same for previous solutions represented by ".valP", ".gradP", and ".hessP" as well as for the test functions "w1" and "w2" and produce arrays, u0[] and w[], respectively.
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
 * \dontinclude mechanoChemStab/userFunctions.cc
 * \skip userFunctions.h
 * \until comp[63]
 *
 *
 * 
 */
