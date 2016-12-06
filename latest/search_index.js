var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#GeometricIntegrators.jl-1",
    "page": "Home",
    "title": "GeometricIntegrators.jl",
    "category": "section",
    "text": "Julia library of geometric integrators for ordinary differential equations and differential algebraic equations.(Image: Build Status) (Image: Coverage Status) (Image: codecov)GeometricIntegrators.jl is a library of geometric integrators for ordinary differential equations and differential algebraic equations in Julia. Its main aim is the implementation and verification of novel geometric integrators, especially with respect to long-time stability and conservation of geometric structures. In order to be able to perform simulations with millions or billions of time steps, the design of the library tries to minimize overhead and maximize performance. For example, all data structures are preallocated and reused so that all runtime allocations are eliminated. GeometricIntegrators.jl provides solvers for various families of integrators as well as facilities to derive such integrators of arbitrary order, e.g., via discrete variational principles."
},

{
    "location": "index.html#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "Currently the following methods are supported or planned to be supported,[x] Explicit Runge-Kutta Methods (ERK),\n[x] Explicit Partitioned Runge-Kutta Methods (EPRK),\n[x] Implicit Partitioned Runge-Kutta Methods (IPRK),\n[x] Fully Implicit Runge-Kutta Methods (FIRK),\n[ ] Diagonally Implicit Runge-Kutta Methods (DIRK),\n[ ] Singly Implicit Runge-Kutta Methods (SIRK),\n[ ] Special Additive Runge-Kutta Methods (SARK),\n[ ] Special Partitioned Additive Runge-Kutta Methods (SPARK),\n[ ] Two-Step Runge-Kutta Methods (TSRK),\n[ ] General Linear Methods (GLM),\n[ ] Splitting Methods (SM).The following families of equations are supported or planned to be supported,[x] Systems of ODEs,\n[ ] Systems of DAEs,\n[x] Partitioned ODEs,\n[ ] Partitioned DAEs,\n[x] Implicit ODEs\n[ ] Implicit DAEswhich can be prescribed manually or obtained as[ ] Euler-Lagrange Equations,\n[ ] Hamilton Equations,\n[ ] Symplectic Equations,\n[ ] Poisson Equations,with[ ] Holonomic Constraints,\n[ ] Nonholonomic Constraints,\n[ ] Dirac Constraints.The following families of integrators are supported or planned to be supported,[ ] Gauss-Legendre Runge-Kutta,\n[ ] Galerkin Variational Integrators,\n[ ] Taylor Variational Integrators,\n[ ] Hamilton-Pontryagin-Galerkin Integrators.Available linear solvers are[x] LU decomposition (LAPACK),\n[ ] LU decomposition (native Julia),\n[ ] Krylov,and nonlinear solvers[ ] Fixed-Point Iteration,\n[ ] Fixed-Point Iteration with Aitken's Acceleration,\n[x] Newton's method,\n[ ] Newton's method with line search,\n[x] Quasi-Newton,either with exact Jacobian or with approximate Jacobian obtained via[x] Finite Differences,\n[ ] Automatic Differentiation,and[ ] Jacobian-free Newton-Krylov."
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Pages = [\"tutorial.md\"]"
},

{
    "location": "index.html#Modules-1",
    "page": "Home",
    "title": "Modules",
    "category": "section",
    "text": "Pages = [\"modules/equations.md\",\n         \"modules/integrators.md\",\n         \"modules/solvers_linear.md\",\n         \"modules/solvers_nonlinear.md\",\n         \"modules/tableaus.md\"\n]"
},

{
    "location": "index.html#License-1",
    "page": "Home",
    "title": "License",
    "category": "section",
    "text": "Copyright (c) 2016 Michael Kraus <michael.kraus@ipp.mpg.de>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
},

{
    "location": "tutorial.html#",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial.html#Tutorial-1",
    "page": "Tutorial",
    "title": "Tutorial",
    "category": "section",
    "text": ""
},

{
    "location": "modules/equations.html#",
    "page": "Equations",
    "title": "Equations",
    "category": "page",
    "text": ""
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.DAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.DAE",
    "category": "Type",
    "text": "DAE: Differential Algebraic Equation\n\nDefines a differential algebraic initial value problem\n\nbeginalign*\ndotq (t) = f(t q(t)) + u(t q(t) lambda(t))   q(t_0) = q_0  \n0 = phi (t q(t) lambda(t))   lambda(t_0) = lambda_0 \nendalign*\n\nwith vector field f, projection u, algebraic constraint phi=0, initial conditions q_0 and lambda_0, the dynamical variable q taking values in mathbbR^m and the algebraic variable lambda taking values in mathbbR^n.\n\nFields\n\nm: dimension of dynamical variable q and the vector field f\nn: dimension of algebraic variable lambda and the constraint phi\nf: function computing the vector field\nu: function computing the projection\nϕ: algebraic constraint\nt₀: initial time\nq₀: initial condition for dynamical variable q\nλ₀: initial condition for algebraic variable lambda\n\nThe function f, providing the vector field, takes three arguments, f(t, q, fq), the functions u and ϕ, providing the projection and the algebraic constraint take four arguments, u(t, q, λ, fu) and ϕ(t, q, λ, fϕ), where t is the current time, q and λ are the current solution vectors, and fq, fu and fϕ are the vectors which hold the result of evaluating the vector field f, the projection u and the algebraic constraint phi on t, q and λ.\n\nExample\n\n    function f(t, q, fq)\n        fq[1] = q[1]\n        fq[2] = q[2]\n    end\n\n    function u(t, q, λ, fu)\n        fu[1] = +λ[1]\n        fu[2] = -λ[1]\n    end\n\n    function ϕ(t, q, λ, fϕ)\n        fϕ[1] = q[2] - q[1]\n    end\n\n    t₀ = 0.\n    q₀ = [1., 1.]\n    λ₀ = [0.]\n\n    dae = DAE(f, u, ϕ, t₀, q₀, λ₀)\n\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.IODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.IODE",
    "category": "Type",
    "text": "IODE: Implicit Ordinary Differential Equation\n\nDefines an implicit initial value problem\n\nbeginalign*\ndotq (t) = v(t)  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) v(t))  \np(t_0) = p_0  \np(t) = g(t q(t) v(t))\nendalign*\n\nwith vector fields f and g, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d. This is a special case of a differential algebraic equation with dynamical variables (qp) and algebraic variable v.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields f and g\nf: function computing the vector field\ng: function determining the algebraic constraint\nt₀: initial time\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions f and g must have the interface\n\n    function f(t, q, v, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nand\n\n    function g(t, q, v, g)\n        g[1] = ...\n        g[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, v is the current velocity and f and g are the vectors which hold the result of evaluating the functions f and g on t, q and v.\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.ODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.ODE",
    "category": "Type",
    "text": "ODE: Ordinary Differential Equation\n\nDefines an initial value problem\n\ndotq (t) = f(t q(t))  qquad q(t_0) = q_0 \n\nwith vector field f, initial condition q_0 and the solution q taking values in mathbbR^d.\n\nFields\n\nd: dimension of dynamical variable q and the vector field f\nf: function computing the vector field\nt₀: initial time\nq₀: initial condition\n\nThe function f providing the vector field must have the interface\n\n    function f(t, q, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q is the current solution vector, and f is the vector which holds the result of evaluating the vector field f on t and q.\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.PDAE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PDAE",
    "category": "Type",
    "text": "Partitioned Differential Algebraic Equation\n\n\n\n"
},

{
    "location": "modules/equations.html#GeometricIntegrators.Equations.PODE",
    "page": "Equations",
    "title": "GeometricIntegrators.Equations.PODE",
    "category": "Type",
    "text": "IODE: Partitioned Ordinary Differential Equation\n\nDefines a partitioned initial value problem\n\nbeginalign*\ndotq (t) = v(t q(t) p(t))  \nq(t_0) = q_0  \ndotp (t) = f(t q(t) p(t))  \np(t_0) = p_0 \nendalign*\n\nwith vector fields v and f, initial conditions (q_0 p_0) and the solution (qp) taking values in mathbbR^d times mathbbR^d.\n\nFields\n\nd: dimension of dynamical variables q and p as well as the vector fields v and f\nv: function computing the vector field v\nf: function computing the vector field f\nt₀: initial time\nq₀: initial condition for q\np₀: initial condition for p\n\nThe functions v and f must have the interface\n\n    function v(t, q, p, v)\n        v[1] = ...\n        v[2] = ...\n        ...\n    end\n\nand\n\n    function f(t, q, p, f)\n        f[1] = ...\n        f[2] = ...\n        ...\n    end\n\nwhere t is the current time, q and p are the current solution vectors and v and f are the vectors which hold the result of evaluating the vector fields v and f on t, q and p.\n\n\n\n"
},

{
    "location": "modules/equations.html#Equations-1",
    "page": "Equations",
    "title": "Equations",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Equations]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/integrators.html#",
    "page": "Integrators",
    "title": "Integrators",
    "category": "page",
    "text": ""
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.DAE,GeometricIntegrators.Integrators.TableauSARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for special additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.Equation,GeometricIntegrators.Integrators.Tableau,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Print error for integrators not implemented, yet.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.IODE,GeometricIntegrators.Integrators.TableauIPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for implicit partitioned Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Integrators.TableauDIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for diagonally implicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Integrators.TableauERK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for explicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Integrators.TableauFIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for fully implicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.ODE,GeometricIntegrators.Integrators.TableauSIRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for singly implicit Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.PDAE,GeometricIntegrators.Integrators.TableauSPARK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for special partitioned additive Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Integrator-Tuple{GeometricIntegrators.Equations.PODE,GeometricIntegrators.Integrators.TableauEPRK,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Integrator",
    "category": "Method",
    "text": "Create integrator for explicit partitioned Runge-Kutta tableau.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorDIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorDIRK",
    "category": "Type",
    "text": "Diagonally implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorEPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorEPRK",
    "category": "Type",
    "text": "Explicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorERK",
    "category": "Type",
    "text": "Explicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorFIRK",
    "category": "Type",
    "text": "Fully implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorIPRK",
    "category": "Type",
    "text": "Implicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorSARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSARK",
    "category": "Type",
    "text": "Special Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSIRK",
    "category": "Type",
    "text": "Singly implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.IntegratorSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.IntegratorSPARK",
    "category": "Type",
    "text": "Special Partitioned Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Solution",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Solution",
    "category": "Type",
    "text": "Create solution for implicit ODE.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Solution",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Solution",
    "category": "Type",
    "text": "Print error for solutions of equations not implemented, yet.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Solution",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Solution",
    "category": "Type",
    "text": "Create solution for DAE.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Solution",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Solution",
    "category": "Type",
    "text": "Create solution for partitioned DAE.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Solution",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Solution",
    "category": "Type",
    "text": "Create solution for partitioned ODE.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Solution",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Solution",
    "category": "Type",
    "text": "Create solution for ODE.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.SolutionDAE",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.SolutionDAE",
    "category": "Type",
    "text": "Solution of a differential algebraic equation.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.SolutionODE",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.SolutionODE",
    "category": "Type",
    "text": "SolutionODE: Solution of an ordinary differential equation\n\nContains all fields necessary to store the solution of an ODE.\n\nFields\n\nnd: dimension of the dynamical variable q\nnt: number of time steps to store\nn0: number of initial conditions\nt:  time steps\nx:  solution x[nd, nt+1, n0] with x[:,0,:] the initial conditions\nntime: number of time steps to compute\nnsave: save every nsave'th time step\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.SolutionPDAE",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.SolutionPDAE",
    "category": "Type",
    "text": "Solution of a partitioned differential algebraic equation.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.SolutionPODE",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.SolutionPODE",
    "category": "Type",
    "text": "Solution of a partitioned ordinary differential equation.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.Tableau",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.Tableau",
    "category": "Type",
    "text": "Holds the information for the various methods' tableaus.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauDIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauDIRK",
    "category": "Type",
    "text": "Holds the tableau of a diagonally implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauEPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauEPRK",
    "category": "Type",
    "text": "Holds the tableau of an explicit partitioned Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauERK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauERK",
    "category": "Type",
    "text": "Holds the tableau of an explicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauFIRK",
    "category": "Type",
    "text": "Holds the tableau of a fully implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauGLM",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauGLM",
    "category": "Type",
    "text": "Holds the tableau of a general linear method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauIPRK",
    "category": "Type",
    "text": "Holds the tableau of a implicit partitioned Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauRK",
    "category": "Type",
    "text": "Holds the tableau of a Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauSARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSARK",
    "category": "Type",
    "text": "Holds the tableau of a spezialized additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauSIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSIRK",
    "category": "Type",
    "text": "Holds the tableau of a singly implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauSPARK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauSPARK",
    "category": "Type",
    "text": "Holds the tableau of a spezialized partitioned additive Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.createHDF5",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.createHDF5",
    "category": "Function",
    "text": "Creates HDF5 file and initialises datasets for ODE solution object.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "Function",
    "text": "Integrate given equation with given tableau for ntime time steps and return solution.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "Function",
    "text": "Apply integrator for ntime time steps and return solution.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate",
    "category": "Function",
    "text": "Integrate ODE specified by vector field and initial condition with given tableau for ntime time steps and return solution.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorDIRK,GeometricIntegrators.Integrators.SolutionODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with diagonally implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorDIRK,GeometricIntegrators.Integrators.SolutionPODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned ODE with diagonally implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorEPRK{T},GeometricIntegrators.Integrators.SolutionPODE{T,tType}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned ODE with explicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorERK{DT,TT,FT},GeometricIntegrators.Integrators.SolutionODE{DT,tType}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with explicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorFIRK,GeometricIntegrators.Integrators.SolutionPODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned ODE with fully implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorFIRK{DT,TT,FT,ST,IT},GeometricIntegrators.Integrators.SolutionODE{DT,TT}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with fully implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorIPRK,GeometricIntegrators.Integrators.SolutionPODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with implicit partitioned Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSARK,GeometricIntegrators.Integrators.SolutionDAE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate DAE with Special Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSIRK,GeometricIntegrators.Integrators.SolutionODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate ODE with singly implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSIRK,GeometricIntegrators.Integrators.SolutionPODE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned ODE with singly implicit Runge-Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.integrate!-Tuple{GeometricIntegrators.Integrators.IntegratorSPARK,GeometricIntegrators.Integrators.SolutionPDAE}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.integrate!",
    "category": "Method",
    "text": "Integrate partitioned DAE with Special Additive Runge Kutta integrator.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.readTableauERKFromFile-Tuple{AbstractString,AbstractString}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.readTableauERKFromFile",
    "category": "Method",
    "text": "Read explicit Runge-Kutta tableau from file.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.showTableau-Tuple{GeometricIntegrators.Integrators.TableauRK{T}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.showTableau",
    "category": "Method",
    "text": "Print Runge-Kutta tableau to standard output.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.writeSolutionToHDF5",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.writeSolutionToHDF5",
    "category": "Function",
    "text": "Append solution to HDF5 file.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.writeSolutionToHDF5-Tuple{GeometricIntegrators.Integrators.Solution,AbstractString}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.writeSolutionToHDF5",
    "category": "Method",
    "text": "Creates HDF5 file, writes solution to file, and closes file.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.writeTableauToFile-Tuple{AbstractString,GeometricIntegrators.Integrators.TableauRK{T}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.writeTableauToFile",
    "category": "Method",
    "text": "Write Runge-Kutta tableau to file.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Tuple{Array{DT,1},Array{TT,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersFIRK{DT,TT,FT}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of fully implicit Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Solvers.function_stages!-Tuple{Array{T,1},Array{T,1},GeometricIntegrators.Integrators.NonlinearFunctionParametersIPRK{T}}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Solvers.function_stages!",
    "category": "Method",
    "text": "Compute stages of implicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersFIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersFIRK",
    "category": "Type",
    "text": "Parameters for right-hand side function of fully implicit Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.NonlinearFunctionParametersIPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.NonlinearFunctionParametersIPRK",
    "category": "Type",
    "text": "Parameters for right-hand side function of implicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauIRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauIRK",
    "category": "Type",
    "text": "Holds the tableau of an implicit Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.TableauPRK",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.TableauPRK",
    "category": "Type",
    "text": "Holds the tableau of a partitioned Runge-Kutta method.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.computeStageP!-Tuple{GeometricIntegrators.Integrators.IntegratorEPRK,Int64,Int64,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.computeStageP!",
    "category": "Method",
    "text": "Compute P stages of explicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.computeStageQ!-Tuple{GeometricIntegrators.Integrators.IntegratorEPRK,Int64,Int64,Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.computeStageQ!",
    "category": "Method",
    "text": "Compute Q stages of explicit partitioned Runge-Kutta methods.\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.getTableauERK4-Tuple{}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.getTableauERK4",
    "category": "Method",
    "text": "Tableau for explicit Runge-Kutta method of order four (1/6 rule)\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.getTableauERK438-Tuple{}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.getTableauERK438",
    "category": "Method",
    "text": "Tableau for explicit Runge-Kutta method of order four (3/8 rule)\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.getTableauExplicitEuler-Tuple{}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.getTableauExplicitEuler",
    "category": "Method",
    "text": "Tableau for explicit Euler method\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.getTableauExplicitMidpoint-Tuple{}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.getTableauExplicitMidpoint",
    "category": "Method",
    "text": "Tableau for explicit midpoint method\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.getTableauHeun-Tuple{}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.getTableauHeun",
    "category": "Method",
    "text": "Tableau for Heun's method\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.getTableauKutta-Tuple{}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.getTableauKutta",
    "category": "Method",
    "text": "Tableau for Kutta's method of order three\n\n\n\n"
},

{
    "location": "modules/integrators.html#GeometricIntegrators.Integrators.readTableauRKHeaderFromFile-Tuple{Any}",
    "page": "Integrators",
    "title": "GeometricIntegrators.Integrators.readTableauRKHeaderFromFile",
    "category": "Method",
    "text": "Reads and parses Tableau metadata from file.\n\n\n\n"
},

{
    "location": "modules/integrators.html#Integrators-1",
    "page": "Integrators",
    "title": "Integrators",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Integrators]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solvers_linear.html#",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solvers_linear.html#Linear-Solvers-1",
    "page": "Linear Solvers",
    "title": "Linear Solvers",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solvers]\nPages   = [\"solvers/linear/linear_solvers.jl\",\n           \"solvers/linear/lu_solver_lapack.jl\"]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/solvers_nonlinear.html#",
    "page": "Nonlinear Solvers",
    "title": "Nonlinear Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "modules/solvers_nonlinear.html#Nonlinear-Solvers-1",
    "page": "Nonlinear Solvers",
    "title": "Nonlinear Solvers",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Solvers]\nPages   = [\"solvers/nonlinear/nonlinear_solvers.jl\",\n           \"solvers/nonlinear/jacobian.jl\",\n           \"solvers/nonlinear/abstract_fixed_point_solver.jl\",\n           \"solvers/nonlinear/fixed_point_solver.jl\",\n           \"solvers/nonlinear/abstract_newton_solver.jl\",\n           \"solvers/nonlinear/newton_solver.jl\",\n           \"solvers/nonlinear/quasi_newton_solver.jl\"]\nOrder   = [:constant, :type, :macro, :function]"
},

{
    "location": "modules/tableaus.html#",
    "page": "Tableaus",
    "title": "Tableaus",
    "category": "page",
    "text": ""
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauERK4-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauERK4",
    "category": "Method",
    "text": "Tableau for explicit Runge-Kutta method of order four (1/6 rule)\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauERK438-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauERK438",
    "category": "Method",
    "text": "Tableau for explicit Runge-Kutta method of order four (3/8 rule)\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauExplicitEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauExplicitEuler",
    "category": "Method",
    "text": "Tableau for explicit Euler method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauExplicitMidpoint-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauExplicitMidpoint",
    "category": "Method",
    "text": "Tableau for explicit midpoint method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauGLRK1-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRK1",
    "category": "Method",
    "text": "Gauss-Legendre Runge-Kutta, s=1\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauGLRK2-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRK2",
    "category": "Method",
    "text": "Gauss-Legendre Runge-Kutta, s=2\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauGLRK3-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauGLRK3",
    "category": "Method",
    "text": "Gauss-Legendre Runge-Kutta, s=3\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauHeun-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauHeun",
    "category": "Method",
    "text": "Tableau for Heun's method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauImplicitEuler-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauImplicitEuler",
    "category": "Method",
    "text": "Implicit Euler\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauImplicitMidpoint-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauImplicitMidpoint",
    "category": "Method",
    "text": "Implicit Midpoint\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauKutta-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauKutta",
    "category": "Method",
    "text": "Tableau for Kutta's method of order three\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauSymplecticEulerA-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSymplecticEulerA",
    "category": "Method",
    "text": "Tableau for symplectic Euler-A method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#GeometricIntegrators.Tableaus.getTableauSymplecticEulerB-Tuple{}",
    "page": "Tableaus",
    "title": "GeometricIntegrators.Tableaus.getTableauSymplecticEulerB",
    "category": "Method",
    "text": "Tableau for symplectic Euler-B method\n\n\n\n"
},

{
    "location": "modules/tableaus.html#Tableaus-1",
    "page": "Tableaus",
    "title": "Tableaus",
    "category": "section",
    "text": "Modules = [GeometricIntegrators.Tableaus]\nOrder   = [:constant, :type, :macro, :function]"
},

]}
