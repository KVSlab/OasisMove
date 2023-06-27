__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2013-11-06"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from dolfin import *
from oasismove.solvers.NSfracStep import *
from oasismove.solvers.NSfracStep import __all__
from petsc4py import PETSc as _PETSc

from oasismove.problems import u_dot_n


def setup(u_components, u, v, p, q, u_w, v_w, bcs, les_model, nut_, scalar_components, V, Q, x_, p_, u_,
          velocity_update_solver, assemble_matrix, homogenize, GradFunction, DivFunction, LESsource, nn_model, nunn_,
          NNsource, mesh, **NS_namespace):
    """Preassemble mass and diffusion matrices.

    Set up and prepare all equations to be solved. Called once, before
    going into time loop.

    """
    # Mass matrix
    M = assemble_matrix(inner(u, v) * dx)

    # Stiffness matrix (without viscosity coefficient)
    K = assemble_matrix(inner(grad(u), grad(v)) * dx)

    # Allocate stiffness matrix for LES that changes with time
    KT = None if les_model == "NoModel" and nn_model == "NoModel" else (
        Matrix(M), inner(grad(u), grad(v)))

    # Pressure Laplacian.
    Ap = assemble_matrix(inner(grad(q), grad(p)) * dx, bcs['p'])

    # Mesh velocity Laplacian
    # alpha = 1 / CellVolume(mesh)
    a_mesh = inner(grad(u_w), grad(v_w)) * dx
    A_mesh = Matrix(assemble_matrix(a_mesh))

    # Allocate coefficient matrix (needs reassembling)
    A = Matrix(M)

    # Allocate Function for holding and computing the velocity divergence on Q
    divu = DivFunction(u_, Q, name='divu',
                       method=velocity_update_solver)

    # Allocate a dictionary of Functions for holding and computing pressure gradients
    gradp = {ui: GradFunction(p_, V, i=i, name='dpd' + ('x', 'y', 'z')[i],
                              bcs=homogenize(bcs[ui]),
                              method=velocity_update_solver)
             for i, ui in enumerate(u_components)}

    # Create dictionary to be returned into global NS namespace
    d = dict(A=A, M=M, K=K, Ap=Ap, a_mesh=a_mesh, A_mesh=A_mesh, divu=divu, gradp=gradp)

    # Allocate coefficient matrix and work vectors for scalars. Matrix differs
    # from velocity in boundary conditions only
    if len(scalar_components) > 0:
        d.update(Ta=Matrix(M))
        if len(scalar_components) > 1:
            # For more than one scalar we use the same linear algebra solver for all.
            # For this to work we need some additional tensors. The extra matrix
            # is required since different scalars may have different boundary conditions
            Tb = Matrix(M)
            bb = Vector(x_[scalar_components[0]])
            bx = Vector(x_[scalar_components[0]])
            d.update(Tb=Tb, bb=bb, bx=bx)

    # Setup for solving convection
    u_ab = as_vector([Function(V) for i in range(len(u_components))])
    a_conv = inner(v, dot(u_ab, nabla_grad(u))) * dx
    a_scalar = a_conv
    LT = None if les_model == "NoModel" else LESsource(
        nut_, u_ab, V, name='LTd')

    NT = None if nn_model == "NoModel" else NNsource(
        nunn_, u_ab, V, name='NTd')

    if bcs['p'] == []:
        attach_pressure_nullspace(Ap, x_, Q)

    d.update(u_ab=u_ab, a_conv=a_conv, a_scalar=a_scalar, LT=LT, KT=KT, NT=NT)
    return d


def get_solvers(use_krylov_solvers, krylov_solvers, krylov_solvers_w, scalar_components, velocity_krylov_solver,
                pressure_krylov_solver, scalar_krylov_solver, mesh_velocity_krylov_solver, **NS_namespace):
    """Return linear solvers.

    We are solving for
       - mesh velocity
       - tentative velocity
       - pressure correction

       and possibly:
       - scalars

    """
    if use_krylov_solvers:
        # tentative velocity solver
        u_prec = PETScPreconditioner(velocity_krylov_solver['preconditioner_type'])
        u_sol = PETScKrylovSolver(velocity_krylov_solver['solver_type'], u_prec)
        u_sol.parameters.update(krylov_solvers)

        # pressure solver
        p_prec = PETScPreconditioner(pressure_krylov_solver['preconditioner_type'])
        p_sol = PETScKrylovSolver(pressure_krylov_solver['solver_type'], p_prec)
        p_sol.parameters.update(krylov_solvers)
        p_sol.set_reuse_preconditioner(True)

        # mesh equation solver
        w_prec = PETScPreconditioner(mesh_velocity_krylov_solver['preconditioner_type'])
        w_sol = PETScKrylovSolver(mesh_velocity_krylov_solver['solver_type'], w_prec)
        w_sol.parameters.update(krylov_solvers_w)

        sols = [u_sol, p_sol, w_sol]
        # scalar solver
        if len(scalar_components) > 0:
            c_prec = PETScPreconditioner(scalar_krylov_solver['preconditioner_type'])
            c_sol = PETScKrylovSolver(scalar_krylov_solver['solver_type'], c_prec)
            c_sol.parameters.update(krylov_solvers)
            sols.append(c_sol)
        else:
            sols.append(None)
    else:
        # tentative velocity solver
        u_sol = LUSolver()
        # pressure solver
        p_sol = LUSolver()
        sols = [u_sol, p_sol]
        # scalar solver
        if len(scalar_components) > 0:
            c_sol = LUSolver()
            sols.append(c_sol)
        else:
            sols.append(None)

    return sols


def assemble_first_inner_iter(A, a_conv, dt, M, scalar_components, les_model, nn_model, a_scalar, K, nu, nut_,
                              nunn_, u_components, LT, KT, NT, b_tmp, b0, x_1, x_2, u_ab, bcs, mesh, boundary, u, v,
                              backflow_facets, backflow_beta, wx_, bw_tmp, bw0, **NS_namespace):
    """Called on first inner iteration of velocity/pressure system.

    Assemble convection matrix, compute rhs of tentative velocity and
    reset coefficient matrix for solve.

    """
    Timer("Assemble first inner iter")
    # Update u_ab used as convecting velocity
    for i, ui in enumerate(u_components):
        u_ab[i].vector().zero()
        u_ab[i].vector().axpy(1.5, x_1[ui])
        u_ab[i].vector().axpy(-0.5, x_2[ui])
        u_ab[i].vector().axpy(-1, wx_[ui])

    A = assemble(a_conv, tensor=A)
    A *= -0.5  # Negative convection on the rhs
    A.axpy(1. / dt, M, True)  # Add mass

    # Set up scalar matrix for rhs using the same convection as velocity
    if len(scalar_components) > 0:
        Ta = NS_namespace['Ta']
        if a_scalar is a_conv:
            Ta.zero()
            Ta.axpy(1., A, True)

    # Add diffusion and compute rhs for all velocity components
    A.axpy(-0.5 * nu, K, True)
    if les_model != "NoModel":
        assemble(nut_ * KT[1] * dx, tensor=KT[0])
        A.axpy(-0.5, KT[0], True)

    if nn_model != "NoModel":
        assemble(nunn_ * KT[1] * dx, tensor=KT[0])
        A.axpy(-0.5, KT[0], True)

    for i, ui in enumerate(u_components):
        # Start with body force
        b_tmp[ui].zero()
        b_tmp[ui].axpy(1., b0[ui])
        # Add transient, convection and diffusion
        b_tmp[ui].axpy(1., A * x_1[ui])
        if les_model != "NoModel":
            LT.assemble_rhs(i)
            b_tmp[ui].axpy(1., LT.vector())
        if nn_model != "NoModel":
            NT.assemble_rhs(i)
            b_tmp[ui].axpy(1., NT.vector())

        bw_tmp[ui].zero()
        bw_tmp[ui].axpy(1., bw0[ui])
    # Reset matrix for lhs
    A *= -1.
    A.axpy(2. / dt, M, True)

    # Add backflow stabilization
    for ID in backflow_facets:
        ds = Measure("ds", domain=mesh, subdomain_data=boundary)
        n = FacetNormal(mesh)
        B = assemble(u_dot_n(u_ab, n) * dot(u, v) * ds(ID))
        as_backend_type(A).mat().axpy(-backflow_beta * 0.5, as_backend_type(B).mat(),
                                      _PETSc.Mat.Structure.SUBSET_NONZERO_PATTERN)
        for i, ui in enumerate(u_components):
            b_tmp[ui].axpy(backflow_beta * 0.5, B * x_1[ui])

    [bc.apply(A) for bc in bcs['u0']]


def attach_pressure_nullspace(Ap, x_, Q):
    """Create null space basis object and attach to Krylov solver."""
    null_vec = Vector(x_['p'])
    Q.dofmap().set(null_vec, 1.0)
    null_vec *= 1.0 / null_vec.norm('l2')
    Aa = as_backend_type(Ap)
    null_space = VectorSpaceBasis([null_vec])
    Aa.set_nullspace(null_space)
    Aa.null_space = null_space


def velocity_tentative_assemble(ui, b, b_tmp, p_, gradp, **NS_namespace):
    """Add pressure gradient to rhs of tentative velocity system."""
    b[ui].zero()
    b[ui].axpy(1., b_tmp[ui])
    gradp[ui].assemble_rhs(p_)
    b[ui].axpy(-1., gradp[ui].rhs)


def velocity_tentative_solve(ui, A, bcs, x_, x_2, u_sol, b, udiff,
                             use_krylov_solvers, **NS_namespace):
    """Linear algebra solve of tentative velocity component."""
    [bc.apply(b[ui]) for bc in bcs[ui]]
    # x_2 only used on inner_iter 1, so use here as work vector
    x_2[ui].zero()
    x_2[ui].axpy(1., x_[ui])
    t1 = Timer("Tentative Linear Algebra Solve")
    u_sol.solve(A, x_[ui], b[ui])
    t1.stop()
    udiff[0] += norm(x_2[ui] - x_[ui])


def pressure_assemble(b, x_, dt, Ap, divu, **NS_namespace):
    """Assemble rhs of pressure equation."""
    divu.assemble_rhs()  # Computes div(u_)*q*dx
    b['p'][:] = divu.rhs
    b['p'] *= (-1. / dt)
    b['p'].axpy(1., Ap * x_['p'])


def pressure_solve(dp_, x_, Ap, b, p_sol, bcs, **NS_namespace):
    """Solve pressure equation."""
    [bc.apply(b['p']) for bc in bcs['p']]
    dp_.vector().zero()
    dp_.vector().axpy(1., x_['p'])
    # KrylovSolvers use nullspace for normalization of pressure
    if hasattr(Ap, 'null_space'):
        p_sol.null_space.orthogonalize(b['p'])

    t1 = Timer("Pressure Linear Algebra Solve")
    p_sol.solve(Ap, x_['p'], b['p'])
    t1.stop()
    # LUSolver use normalize directly for normalization of pressure
    if bcs['p'] == []:
        normalize(x_['p'])

    dpv = dp_.vector()
    dpv.axpy(-1., x_['p'])
    dpv *= -1.


def velocity_update(u_components, bcs, gradp, dp_, dt, x_, **NS_namespace):
    """Update the velocity after regular pressure velocity iterations."""
    for ui in u_components:
        gradp[ui](dp_)
        x_[ui].axpy(-dt, gradp[ui].vector())
        [bc.apply(x_[ui]) for bc in bcs[ui]]


def mesh_velocity_assemble(A_mesh, ui, bw, bw_tmp, a_mesh, bc_mesh, A_cache, **NS_namespace):
    """Assemble rhs of mesh velocity equation."""
    A_mesh.zero()
    A_mesh.axpy(1, A_cache[(a_mesh, tuple(bc_mesh[ui]))], True)
    bw[ui].zero()
    bw[ui].axpy(1., bw_tmp[ui])


def mesh_velocity_solve(A_mesh, bw, wx_, w_, dof_map, dt, coordinates, w_sol, ui, bc_mesh, OasisTimer, **NS_namespace):
    """Solve mesh equation."""
    [bc.apply(bw[ui]) for bc in bc_mesh[ui]]
    w_sol.solve(A_mesh, wx_[ui], bw[ui])

    arr = w_[ui].vector().get_local(dof_map)
    mesh_tolerance = 1e-15
    if mesh_tolerance < abs(arr.min()) + abs(arr.max()):
        coordinates[:, int(ui[-1])] += arr * dt


def scalar_assemble(a_scalar, a_conv, Ta, dt, M, scalar_components, Schmidt_T, KT, nu, Schmidt, b, K, x_1, b0,
                    les_model, nn_model, **NS_namespace):
    """Assemble scalar equation."""
    # Just in case you want to use a different scalar convection
    if a_scalar is not a_conv:
        assemble(a_scalar, tensor=Ta)
        Ta *= -0.5  # Negative convection on the rhs
        Ta.axpy(1. / dt, M, True)  # Add mass

    # Compute rhs for all scalars
    for ci in scalar_components:
        # Add diffusion
        Ta.axpy(-0.5 * nu / Schmidt[ci], K, True)
        if les_model != "NoModel":
            Ta.axpy(-0.5 / Schmidt_T[ci], KT[0], True)
        if nn_model != "NoModel":
            Ta.axpy(-0.5 / Schmidt_T[ci], KT[0], True)

        # Compute rhs
        b[ci].zero()
        b[ci].axpy(1., Ta * x_1[ci])
        b[ci].axpy(1., b0[ci])

        # Subtract diffusion
        Ta.axpy(0.5 * nu / Schmidt[ci], K, True)
        if les_model != "NoModel":
            Ta.axpy(0.5 / Schmidt_T[ci], KT[0], True)
        if nn_model != "NoModel":
            Ta.axpy(0.5 / Schmidt_T[ci], KT[0], True)

    # Reset matrix for lhs - Note scalar matrix does not contain diffusion
    Ta *= -1.
    Ta.axpy(2. / dt, M, True)


def scalar_solve(ci, scalar_components, Ta, b, x_, bcs, c_sol, nu, Schmidt, K, **NS_namespace):
    """Solve scalar equation."""

    Ta.axpy(0.5 * nu / Schmidt[ci], K, True)  # Add diffusion
    if len(scalar_components) > 1:
        # Reuse solver for all scalars. This requires the same matrix and vectors to be used by c_sol.
        Tb, bb, bx = NS_namespace['Tb'], NS_namespace['bb'], NS_namespace['bx']
        Tb.zero()
        Tb.axpy(1., Ta, True)
        bb.zero()
        bb.axpy(1., b[ci])
        bx.zero()
        bx.axpy(1., x_[ci])
        [bc.apply(Tb, bb) for bc in bcs[ci]]
        c_sol.solve(Tb, bx, bb)
        x_[ci].zero()
        x_[ci].axpy(1., bx)

    else:
        [bc.apply(Ta, b[ci]) for bc in bcs[ci]]
        c_sol.solve(Ta, x_[ci], b[ci])
    Ta.axpy(-0.5 * nu / Schmidt[ci], K, True)  # Subtract diffusion
