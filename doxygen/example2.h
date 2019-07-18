/**
 * @page example2 Example 2 : Cahn-Hilliard (one species)
 * \dontinclude CahnHilliard_oneSpecies/2D/userFunctions.cc
 *
 * This example implements the Cahn-Hilliard equation for phase-field modeling of a single species,
 * as described by the following weak form of the PDE. The scalar field is composition, \f$c\f$.
 * Note the application of the higher-order Dirichlet boundary condition \f$\nabla c\cdot\boldsymbol{n}=0\f$
 * using Nitsche's method.
 *
 * Cahn-Hilliard:
 *
 * \f{eqnarray*}{
 * 0 &=& \int_\Omega \left(w_1\frac{c - c_{prev}}{\mathrm{d}t} + 
 * M\left(\nabla w_1\cdot(f_{,cc}\nabla c) + \kappa_1\nabla^2 w_1\nabla^2 c\right)\right) dV\\
 * &\phantom{=}& - \int_{\partial\Omega} \left(w_1j_n + 
 * M\kappa_1\left(\nabla^2c(\nabla w_1\cdot\boldsymbol{n}) + \nabla^2w_1(\nabla c\cdot\boldsymbol{n})\right)
 *  - \tau(\nabla w_1\cdot\boldsymbol{n})(\nabla c\cdot\boldsymbol{n})\right) dS
 * \f}
 *
 * Free energy density:
 *
 * \f{eqnarray*}{
 * f(c) = \alpha(c - c_a)^2(c - c_b)^2
 * \f}
 *
 * The current settings prescribe random initial conditions and zero-flux
 * boundary conditions. With these settings, the following evolution 
 * of the concentration is obtained:
 *
 * \htmlonly <style>div.image img[src="example2.png"]{width:20cm;}</style> \endhtmlonly
 * @image html example2.png 
 *
 * Implementation: Level 1 users
 * ==============================
 *
 * To implement this model, we will specify the following through defining user functions: <br>
 * - Initial conditions <br>
 * - Constitutive model (via free energy density functions) <br>
 * - Parameter values <br>
 * - Weak form of the PDEs <br>
 *
 * First, we include the header file declaring the required user functions. These functions will be defined in this file.
 *
 * \line userFunctions
 *
 * Now, we first define any optional user functions. Optional user functions have a default definition that can be 
 * redefined by the user using a function pointer.
 * This will be done in the \c defineParameters function. The available list of optional user functions includes:
 * \c boundaryConditions, \c scalarInitialConditions, \c vectorInitialConditions, \c loadStep, \c adaptiveTimeStep, and \c projectFields.
 * In this example, we redefine only the \c scalarInitialConditions function, while using the default functions for the others.
 *
 * <b> The \c scalarInitialConditions function </b>
 *
 * We initialized the composition field to be random about 0.5.
 *
 * \skip template
 * \until //end
 *
 * <b> Free energy density derivative functions </b>
 *
 * This phase-field implementation requires the second derivative of the chemical free energy density function
 * \f$f(c,\eta) = \alpha(c - c_a)^2(c - c_b)^2\f$. We define the function computing
 * \f$\partial^2 f/\partial c^2\f$ here.
 * Note that this free energy derivative function is used only in this file. It is not a member of any class,
 * nor will we use it to set any function pointers.
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
 * Many of these values can be overwritten by the parameters.prm file, which we will look at later.
 *
 * \skip template
 * \until void
 *
 * Here, we define the mesh by setting the number of elements in each direction, e.g. a 100x100 element mesh.
 *
 * \skip user.N[0]
 * \until user.N[1]
 *
 * We also define the dimensions of the domain, e.g. a unit square.
 *
 * \skip user.L[0]
 * \until user.L[1]
 *
 * We can define a periodic (or partially periodic) domain. The default is no periodicity in all directions.
 * Here, we override the default and define periodicity in the x direction.
 *
 * \skip user.periodic[0]
 * \until user.periodic[0]
 *
 * We can define additional material parameters that are not explicity listed in the \c user structure by
 * defining elements of the \c matParam C++ map, which maps \c std::string to \c double. These values can also be overwritten
 * in the parameters file.
 *
 * \skip "inFlux"
 * \until "c_b"
 *
 * We define the initial time step and total simulation time. We also have the options to use restart files, in which case
 * we would set the iteration index and time at which to start. We leave these values at zero to begin a new simulation.
 * We also have the option to output results at regular intervals (e.g. every 5 time steps).
 *
 * \skip user.dtVal
 * \until user.skipOutput
 *
 * We specify the number of vector and scalar solution and projection fields by adding the name of each field to
 * their respective vector. Here, we have one scalar solution field (the composition).
 * We do not use any vector solution fields or projection fields in this example.
 *
 * \line "c"
 * 
 * We can specify the polynomial order of the basis splines, as well as the global continuity.
 * Note that the global continuity must be less than the polynomial order.
 * Here, we use quadratic basis functions with C-1 global continuity.
 *
 * \skip polyOrder
 * \until globalContinuity
 *
 * Finally, we redirect the desired user function pointers to the \c scalarInitialConditions function that we
 * defined above. This completes the \c defineParameters function.
 *
 * \skip scalarInitialConditions
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
 * \c trans( ) - transpose 2nd order tensor \n
 * \c trace( ) - trace of 2nd order tensor \n
 * \c det( ) - determinant of 2nd order tensor \n
 * \c inv( ) - inverse of 2nd order tensor \n
 *
 *
 * The example code here implements the weak form for the Cahn-Hilliard equation, as shown above.
 *
 * First, we set the values for necessary parameters, using some predefined material parameters.
 *
 * \skip dt
 * \until tau
 *
 * Next, we get the values for the free energy derivative \f$f_{,cc}\f$  based on the current quadrature point.
 *
 * \skip f_cc;
 * \until f_cc =
 *
 * Now, we compute the residual in a manner very similar to the analytical form:
 *
 * \f{eqnarray*}{
 * 0 &=& \int_\Omega \left(w_1\frac{c - c_{prev}}{\mathrm{d}t} + 
 * M\left(\nabla w_1\cdot(f_{,cc}\nabla c) + \kappa_1\nabla^2 w_1\nabla^2 c\right)\right) dV\\
 * &\phantom{=}& - \int_{\partial\Omega} \left(w_1j_n + 
 * M\kappa_1\left(\nabla^2c(\nabla w_1\cdot\boldsymbol{n}) + \nabla^2w_1(\nabla c\cdot\boldsymbol{n})\right)
 *  - \tau(\nabla w_1\cdot\boldsymbol{n})(\nabla c\cdot\boldsymbol{n})\right) dS
 * \f}
 *
 * \skip r =
 * \until //end
 *
 * Finally, we include a file that instatiates the template functions \c defineParameters and \c residual. This bit of code
 * will generally be the same for any problem (unless you decide to use a different automatic differentation library);
 * the user does not need to modify it.
 *
 * \line "userFunctionsInstantiation.h"
 *
 * The complete implementation can be found at  <a href="https://github.com/mechanoChem/mechanoChemIGA/blob/master/initBounValProbs/CahnHilliard_oneSpecies/2D/userFunctions.cc">Github</a>.
 *
 * Parameters file: Interface for level 2 users
 * ==============================
 *
 * Now let's look at the parameters file, \c parameters.prm. The advantages of the parameters file are that
 * these values can be changed without recompiling the code and it can provide a clean interface to the code.
 * \dontinclude CahnHilliard_oneSpecies/2D/parameters.prm
 *
 * The parameters defined in the parameters file overwrite any previous values defined in the \c defineParameters function.
 * Anything following the pound sign (#) is a comment. A parameter is defined using the syntax: 
 *
 * \c set \c parameterName \c = \c parameterValue
 *
 * There is a set list of variables that can be read from the parameters file. Anything else will be added to 
 * the \c matParam structure with a double number type. Tensor objects can follow the format: 1 x 1 or [1,1] or (1,1), 
 * where the number of components must equal the spatial dimension of the problem.
 *
 * In this example file, we begin by specifying the spatial dimension, the geometry dimensions, and the mesh size:
 *
 * \skip dim
 * \until set N
 *
 * Next, we define some parameters that are specific to this problem,
 * so they become elements of \c matParam (see the \c residual and \defineParameters functions above).
 *
 * \skip Free energy
 * \until kappa
 *
 * We then define time stepping, restart information, output frequency, and spline parameters.
 *
 * \skip Time stepping
 * \until globalContinuity
 *
 * Note that we don't need to include all (or even any) of these parameters in this file. We defined default values previously. 
 *
 * The complete parameters file can be found at  <a href="https://github.com/mechanoChem/mechanoChemIGA/blob/master/initBounValProbs/CahnHilliard_oneSpecies/2D/parameters.prm">Github</a>.
 *
 */
