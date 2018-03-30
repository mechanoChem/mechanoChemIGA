/**
 * @page example1 Example 1 : Nongradient, finite strain mechanics
 * \dontinclude nonGradientMechanics/3D/userFunctions.cc
 *
 * To model nongradient, finite strain elasticity, we will specify the following through defining user functions: <br>
 * - Boundary conditions <br>
 * - Derived fields for output (e.g. eqivalent stress) <br>
 * - Constitutive model (via the 1st Piola-Kirchhoff stress) <br>
 * - Parameter values <br>
 * - Weak form of the PDE <br>
 *
 * First, we include the header file declaring the required user functions. These functions will be defined in this file.
 *
 * \line userFunctions
 *
 * Now, we first define any optional user functions. Optional user functions have a default definition that can be 
 * redefined by the user using a function pointer.
 * This will be done in the \c defineParameters function. The available list of optional user functions includes:
 * \c boundaryConditions, \c scalarInitialConditions, \c vectorInitialConditions, \c loadStep, \c adaptiveTimeStep, and \c projectFields.
 * In this example, we redefine the \c boundaryConditions and \c projectFields functions, while using the default functions for the others.
 *
 * <b> The \c boundaryConditions function </b>
 *
 * This function defines Dirichlet boundary conditions using PetIGA's \c IGASetBoundaryValue function.
 * The arguments to this function are as follows: the iga object (user.iga),
 * the "axis" (0, 1, or 2, corresponding to the x, y, or z-axis),
 * the "side" (0 or 1), the "dof", and the "value" that is to be imposed.
 * Note that this can only set a uniform value for a degree-of-freedom on any side.
 * Here, we fix all three degrees-of-freedom on the surface at z=0.
 *
 * \skip template
 * \until //end
 *
 * <b> The \c projectFields function </b>
 *
 * If there are field values derived from the solution fields that are of interest, we can compute these
 * values at each quadrature point and project the value to the nodes. Here, we compute the equivalent stress 
 * using the 1st Piola-Kirchhoff stress. Scalar values are stored in the \c scalarProjections vector and
 * vector values are stored in the \c vectorProjections vector.
 *
 * \skip template
 * \until //end
 *
 * <b> The \c get1stPiolaKirchhoff function </b>
 *
 * This function defines the 1st Piola-Kirchhoff stress. However, it is used only in this file 
 * (by the \c projectFields and \c residual functions), so it is not a class member function nor does it have an associated function pointer.
 * Additional nonmember functions can be defined in this file, if they are only used by other functions in this same file.
 *
 * \skip template
 * \until //end
 *
 * <b> The \c defineParameters function </b>
 *
 * The user is required to define the \c defineParameters and \c residual functions. The \c defineParameters defines variables
 * and functions in the \c AppCtx object. The \c AppCtx object is defined
 * in the appCtx.h file. This function is used to define any values in \c user that will be needed in the problem.
 * It is also used to set any function pointers for user functions that we have redefined.
 *
 * \skip template
 * \until void
 *
 * Here, we define the mesh by setting the number of elements in each direction, e.g. a 10x10x10 element mesh.
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
 * their respective vector. Here, we have one vector solution field (the displacement) and one scalar projection field
 * (the von Mises stress). We do not use any scalar solution or vector projection fields in this example.
 *
 * \skip "displacement"
 * \until "vonMises"
 * 
 * We can specify the polynomial order of the basis splines, as well as the global continuity.
 * Note that the global continuity must be less than the polynomial order.
 * Here, we use a linear basis function with C-0 continuity.
 *
 * \skip polyOrder
 * \until globalContinuity
 *
 * We now define the 4th order elasticity tensor.
 *
 * \skip double E
 * \until //end
 *
 * Finally, we redirect the desired user function pointers to the \c boundaryConditions and \c projectFields functions that we
 * defined above. This completes the \c defineParameters function.
 *
 * \skip boundaryConditions
 * \until //end
 *
 * <b> The \c residual function </b>
 *
 * The residual function defines the residual that is to be driven to zero.
 * This is the central function of the code.
 * It is set up to follow the analytical weak form of the PDE.
 * It has a number of arguments that give problem information at the current quadrature point.
 *
 * \skip template
 * \until &r
 *
 * \c dV is a boolean, "true" if \c residual is being called for the volume integral and 
 * "false" if \c residual is being called for the surface integral.\n
 * \c dS is a boolean, "false" if \c residual is being called for the volume integral and 
 * "true" if \c residual is being called for the surface integral.\n
 * \c x gives the coordinates of the quadrature point.\n
 * \c normal gives the unit normal for a surface quadrature point.\n
 * \c c gives the information (values, gradients, etc.) for the scalar solution fields at the current quadrature point
 * (see documentation for solutionScalars class).\n
 * \c u gives the information (values, gradients, etc.) for the vector solution fields at the current quadrature point
 * (see documentation for solutionVectors class).\n
 * \c w1 gives the information for the scalar test functions.\n
 * \c w2 gives the information for the vector test functions.\n
 * \c user is a structure available for parameters related to the initial boundary value problem (e.g. elasticity tensor).\n
 * \c r stores the scalar value of the residual for the weak form of the PDE which is then used by the core assembly functions.
 *
 * The following functions are available for the solution objects \c c and \c u,
 * where the argument is the field index, i.
 *
 * \c c.val(i) - Value of scalar field i, scalar \n
 * \c c.grad(i) - Gradient of scalar field i, 1st order tensor \n
 * \c c.hess(i) - Hessian of scalar field i, 2nd order tensor \n
 * \c c.laplacian(i) - Laplacian of scalar field i, scalar \n
 * \c c.valP(i) - Value of scalar field i at previous time step, scalar \n
 * \c c.gradP(i) - Gradient of scalar field i at previous time step, 1st order tensor \n
 * \c c.hessP(i) - Hessian of scalar field i at previous time step, 2nd order tensor \n
 * \c c.laplacianP(i) - Laplacian of scalar field i at previous time step, scalar
 *
 * \c u.val(i) - Value of vector field i, 1st order tensor \n
 * \c u.grad(i) - Gradient of vector field i, 2nd order tensor \n
 * \c u.hess(i) - Hessian of vector field i, 3rd order tensor \n
 * \c u.valP(i) - Value of vector field i at previous time step, 1st order tensor \n
 * \c u.gradP(i) - Gradient of vector field i at previous time step, 2nd order tensor \n
 * \c u.hessP(i) - Hessian of vector field i at previous time step, 3rd order tensor
 *
 * Similar functions are available for the test functions. Also, the following tensor operations are useful:
 *
 * Tensor operations: \n
 * \c operator+ - tensor addition \n
 * \c operator- - tensor subraction \n
 * \c operator* - single contraction between tensors or scalar multiplication \n
 * \c double_contract - double contraction of two 2nd order tensors, or
 *                   a 4th order tensor and a 2nd order tensor. \n
 * \c trans( ) - transpose 2nd order tensor
 * \c trace( ) - trace of 2nd order tensor
 * \c det( ) - determinant of 2nd order tensor
 * \c inv( ) - inverse of 2nd order tensor
 *
 * The example code here implements the weak form for finite strain elasticity,
 * \f$\int_\Omega (\nabla{\boldsymbol{w}}:\boldsymbol{P}) dV - \int_{\partial\Omega} (\boldsymbol{w}\cdot\boldsymbol{h}) dS = 0\f$
 * with the Neumann boundary condtion \f$\boldsymbol{h} = \langle 0,0,1.e9x_1\rangle\f$ on \f$x_3=1\f$.
 * First, we get the values for \f$\boldsymbol{P}\f$ and \f$\boldsymbol{h}\f$, based on the current quadrature point.
 *
 * \skip //Elasticity
 * \until h[2]
 *
 * Now, we compute the residual in a manner very similar to the analytical form
 *
 * \skip \int_\Omega
 * \until //end
 *
 * Finally, we include a file that instatiates the template functions \c defineParameters and \c residual. This bit of code
 * will generally be the same for any problem (unless you decide to use a different automatic differentation library),
 * the user does not need to modify it.
 *
 * \line "userFunctionsInstantiation.h"
 *
 * The complete code
 * ==============================
 *
 * \include nonGradientMechanics/3D/userFunctions.cc
 */
