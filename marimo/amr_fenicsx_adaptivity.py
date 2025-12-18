# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "marimo==0.18.4",
#   "numpy",
#   "matplotlib",
#   "mpi4py",
#   "petsc4py",
#   "fenics-dolfinx==0.10.0",
#   "trame-vtk",  # needed for PyVista HTML export
#   "trimesh",    # fallback inline renderer if trame-vtk is unavailable
# ]
# ///
"""
Adaptive mesh refinement (AMR) demo in FEniCSx (DOLFINx) + marimo.

Run:
  marimo edit amr_fenicsx_adaptivity.py
or:
  marimo run  amr_fenicsx_adaptivity.py

Notes
- This tutorial implements a simple *residual-based* indicator and a
  mark/refine loop. It is meant for teaching, not for production verification.
- Works in serial or MPI. Plots only on rank 0.
"""

import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    # Adaptive mesh refinement (AMR) in FEniCSx (DOLFINx)

    We solve a Poisson problem on the unit square with a **localized Gaussian source**
    so the solution develops a steep feature near the source. We then:

    1. solve on a mesh,
    2. build a **residual-based local indicator** on each cell,
    3. **mark** a fraction of cells with the largest indicators,
    4. **refine** the mesh and repeat.

    This mimics the classical adaptive loop: **SOLVE → ESTIMATE → MARK → REFINE**.
    """)
    return


@app.cell
def _():
    import os
    from pathlib import Path

    import numpy as np
    import matplotlib.pyplot as plt

    from mpi4py import MPI
    from petsc4py import PETSc

    # FEniCSx JIT caches generated code (FFCx) under XDG_CACHE_HOME. Always use
    # a local `.cache/` next to this notebook to avoid permission issues.
    local_cache = Path(__file__).resolve().parent / ".cache"
    local_cache.mkdir(parents=True, exist_ok=True)
    os.environ["XDG_CACHE_HOME"] = str(local_cache)

    import ufl
    from dolfinx import fem, mesh, plot as dplot
    from dolfinx.fem import petsc as fem_petsc
    return MPI, PETSc, dplot, fem, fem_petsc, mesh, np, plt, ufl


@app.cell
def _(MPI):
    comm = MPI.COMM_WORLD
    rank = comm.rank
    return comm, rank


@app.cell
def _(mo):
    mo.md(r"""
    ## Problem setup

    We solve

    \[
    -\nabla\cdot\bigl(\kappa \nabla u\bigr) = f \quad \text{in } \Omega=(-1,1)^2, \qquad
    u=0 \quad \text{on } \partial \Omega.
    \]

    with a discontinuous diffusion coefficient \(\kappa(x)\) and a constant right-hand side \(f=1\).
    """)
    return


@app.cell
def _(PETSc, fem, mesh, np, ufl):
    def build_problem(domain: mesh.Mesh, degree: int = 1):
        """Return (V, kappa_expr, kappa, f_rhs, a, L, bcs) for diffusion with homogeneous Dirichlet BC."""
        V = fem.functionspace(domain, ("Lagrange", degree))

        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        x = ufl.SpatialCoordinate(domain)
        # Discontinuous indicator defining an inhomogeneous diffusion coefficient:
        #   chi(x,y) = 1 if (x*y > 0) AND (|x| < 0.5) AND (|y| < 0.5), else 0
        mask_xy = ufl.conditional(ufl.gt(x[0] * x[1], 0.0), 1.0, 0.0)
        mask_x = ufl.conditional(ufl.lt(abs(x[0]), 0.5), 1.0, 0.0)
        mask_y = ufl.conditional(ufl.lt(abs(x[1]), 0.5), 1.0, 0.0)
        kappa_expr = mask_xy * mask_x * mask_y

        # "Heat/diffusion coefficient"
        kappa = 1.0 + 20.0 * kappa_expr

        # Forcing active only where kappa_expr == 1 (same support as jump in kappa)
        f_rhs = kappa_expr * fem.Constant(domain, PETSc.ScalarType(1.0))

        a = ufl.inner(kappa * ufl.grad(u), ufl.grad(v)) * ufl.dx
        L = f_rhs * v * ufl.dx

        # Homogeneous Dirichlet BC on the whole boundary.
        def boundary(x):
            # "near" is not required; for unit square with float coords, this is ok.
            return (
                np.isclose(x[0], -1.0)
                | np.isclose(x[0], 1.0)
                | np.isclose(x[1], -1.0)
                | np.isclose(x[1], 1.0)
            )

        facets = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, boundary)
        dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, facets)
        bc = fem.dirichletbc(PETSc.ScalarType(0), dofs, V)

        return V, kappa_expr, kappa, f_rhs, a, L, [bc]
    return (build_problem,)


@app.cell
def _(comm, mesh, np):
    # Start with a coarse mesh on (-1, 1)^2.
    domain0 = mesh.create_rectangle(
        comm,
        [np.array([-1.0, -1.0]), np.array([1.0, 1.0])],
        [20, 20],
        cell_type=mesh.CellType.triangle,
    )
    return (domain0,)


@app.cell
def _(build_problem, fem_petsc):
    def solve_poisson(domain):
        V, kappa_expr, kappa, f_rhs, a, L, bcs = build_problem(domain, degree=1)

        problem = fem_petsc.LinearProblem(
            a, L, bcs=bcs,
            petsc_options_prefix="poisson_",
            petsc_options={
                "ksp_type": "cg",
                "pc_type": "hypre",
                "ksp_rtol": 1.0e-10,
            },
        )
        uh = problem.solve()
        uh.name = "u"
        return V, kappa_expr, kappa, f_rhs, uh
    return (solve_poisson,)

@app.cell
def _(mo):
    mo.md(
        r"""
## 1) Error equation and residual functional

Continuous weak form (exact solution \(u\in H_0^1(\Omega)\)):

\[
(\nabla u,\nabla v) = (f,v)\quad \forall v\in H_0^1(\Omega).
\]

Subtract the discrete term \((\nabla u_h,\nabla v)\) from both sides and set \(e:=u-u_h\):

\[
(\nabla e,\nabla v) = (f,v) - (\nabla u_h,\nabla v).
\]

Define the residual functional

\[
R(v) := (f,v) - (\nabla u_h,\nabla v),
\]

so we have the compact identity

\[
\boxed{(\nabla e,\nabla v) = R(v)\quad \forall v\in H_0^1(\Omega).}
\]
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 2) Galerkin orthogonality (compact)

For all \(v_h\in V_h\subset H_0^1(\Omega)\),

\[
(\nabla u,\nabla v_h)=(f,v_h),\qquad (\nabla u_h,\nabla v_h)=(f,v_h).
\]

Subtracting gives

\[
\boxed{(\nabla e,\nabla v_h)=0\quad \forall v_h\in V_h}
\qquad\Longleftrightarrow\qquad
\boxed{R(v_h)=0\quad \forall v_h\in V_h.}
\]

This will let us replace \(v\) by \(v-I_hv\) later without changing \(R(\cdot)\).
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 3) Elementwise residual form and the \(f-f_h\) replacement

Start from

\[
R(v)=(f,v)-(\nabla u_h,\nabla v)
=
\sum_{T\in\mathcal T_h}
\Big[(f,v)_T-(\nabla u_h,\nabla v)_T\Big].
\]

Integrate by parts on each element \(T\):

\[
(\nabla u_h,\nabla v)_T
=
-(\Delta u_h,v)_T
+
\langle \nabla u_h\cdot n_T,\ v\rangle_{\partial T}.
\]

Combining these and joining boundary contributionsion from neighboring elements into jump term:

\[
R(v)
=
\sum_T (f+\Delta u_h,\ v)_T
\;-\;
\sum_T \langle \nabla u_h\cdot n_T,\ v\rangle_{\partial T}
=
\sum_T (R_T,\ v)_T
\;-\;
\sum_{E\in\mathcal E_h^{\rm int}} \langle J_E,\ v\rangle_E,
\]

where

\[
R_T := f+\Delta u_h,
\qquad
J_E := [\nabla u_h\cdot n]_E.
\]

"""
    )
    return

@app.cell
def _(mo):
    mo.md(
        r"""
## 4) Local residual bounds (cellwise viewpoint)

We now estimate the residual **locally on each cell**, without yet invoking the
energy norm of the error.

Using Galerkin orthogonality \(R(v_h)=0\) for all \(v_h\in V_h\),
choose an interpolation \(I_h v\in V_h\) and write

\[
R(v)
=
R(v-I_hv).
\]

Insert the elementwise residual representation:

\[
R(v-I_hv)
=
\sum_T (R_T,\ v-I_hv)_T
\;-\;
\sum_{E\in\mathcal E_h^{\rm int}}
\langle J_E,\ v-I_hv\rangle_E.
\]

Consider first a **single cell contribution**.
By Cauchy–Schwarz,

\[
|(R_T,\ v-I_hv)_T|
\le
\|R_T\|_{L^2(T)}\,\|v-I_hv\|_{L^2(T)}
\lesssim
h_T\,\|R_T\|_{L^2(T)}\,\|\nabla v\|_{L^2(\omega_T)},
\]

where we have applied  a standard interpolation estimate:

\[
\|v-I_hv\|_{L^2(T)}
\lesssim
h_T\,\|\nabla v\|_{L^2(\omega_T)}.
\]


A completely analogous argument using trace inequalities yields, for interior facets,

\[
|\langle J_E,\ v-I_hv\rangle_E|
\lesssim
h_E^{1/2}\,\|J_E\|_{L^2(E)}\,\|\nabla v\|_{L^2(\omega_E)}.
\]

At this stage, everything is local: one cell, one facet, one gradient of \(v\).
"""
    )
    return






@app.cell
def _(mo):
    mo.md(
        r"""
## 5) From local bounds to the cellwise energy error estimator

Sum the local bounds from Step 4 over all elements and facets, and apply
Cauchy–Schwarz across the mesh:

\[
|R(v-I_hv)|
\;\lesssim\;
\Bigg(
\sum_T h_T^2\|R_T\|_{L^2(T)}^2
+
\sum_{E\in\mathcal E_h^{\rm int}} h_E\|J_E\|_{L^2(E)}^2
\Bigg)^{1/2}
\,
\|\nabla v\|_{L^2(\Omega)}.
\]

Now invoke the dual characterization of the energy norm:

\[
\|\nabla e\|_{L^2(\Omega)}
=
\sup_{v\in H_0^1(\Omega),\,v\neq 0}
\frac{R(v)}{\|\nabla v\|_{L^2(\Omega)}}
=
\sup_{v\neq 0}
\frac{R(v-I_hv)}{\|\nabla v\|_{L^2(\Omega)}}.
\]

Insert the bound above.
The factor \(\|\nabla v\|\) cancels, yielding

\[
\|\nabla e\|_{L^2(\Omega)}^2
\;\lesssim\;
\Bigg(
\sum_T h_T^2\|R_T\|_{L^2(T)}^2
+
\sum_{E\in\mathcal E_h^{\rm int}} h_E\|J_E\|_{L^2(E)}^2
\Bigg).
\]

This motivates the **cellwise residual indicator**

\[
\boxed{
\eta_T^2
:=
h_T^2\|R_T\|_{L^2(T)}^2
+
\frac12
\sum_{E\subset\partial T\cap\mathcal E_h^{\rm int}}
h_E\|J_E\|_{L^2(E)}^2.
}
\]

*Side comment (lecturer): The factor \(1/2\) simply splits each interior facet
between its two neighboring cells — this is what makes marking robust.*

Summing \(\eta_T^2\) over all cells recovers the global estimator and yields
a reliable bound for the energy norm of the error.



===
---
TODO: postpone 
### Data oscillation and computability

If the right-hand side \(f\) is not exactly representable in the discrete setting,
introduce a local polynomial approximation \(f_h\)
(e.g. the elementwise \(L^2\)-projection).

Add and subtract \(f_h\) inside the element residual:

\[
R_T
=
(f-f_h) + (f_h+\Delta u_h).
\]

*Side comment (lecturer): Note carefully — \(f_h+\Delta u_h\) is **not zero in general**.
The discrete problem only enforces \((f_h+\Delta u_h,v_h)_T=0\) for test functions
\(v_h\in V_h\), not pointwise or for arbitrary \(v\).*

Thus the residual splits into

- a **computable discrete residual** \(f_h+\Delta u_h\),
- a **data oscillation term** \(f-f_h\).

This distinction will matter later, but for now we keep the exact form \(R_T\).


"""
    )
    return



@app.cell
def _(fem, fem_petsc, ufl):
    def project_to_space(expr, V, prefix="proj_"):
        """
        L2 projection: find p in V such that (p, v) = (expr, v) for all v in V.
        Used to project continuous data (f) to the discrete space.
        """
        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)
        problem = fem_petsc.LinearProblem(
            ufl.inner(u, v) * ufl.dx,
            ufl.inner(expr, v) * ufl.dx,
            petsc_options_prefix=prefix,
            petsc_options={"ksp_type": "cg", "pc_type": "jacobi", "ksp_rtol": 1.0e-10},
        )
        p = fem.Function(V, name="proj")
        p.x.array[:] = problem.solve().x.array
        return p

    return (project_to_space,)


@app.cell
def _(PETSc, fem, fem_petsc, np, project_to_space, ufl):
    def cell_residual_indicator(domain, kappa, f_rhs, uh):
        """
        Compute a per-cell indicator eta_T using a residual-style estimator:

        - Cell residual: r = f_rhs - f_h where f_h is the L2 projection of f_rhs onto V
        - Interior facet flux jump: [[kappa*grad(uh)·n]]

        We assemble the per-cell quantity into a DG0 vector (one DOF per cell).
        Facet contributions are split between neighboring cells via `avg(w)`.
        """
        # DG0 space: one DOF per cell
        W = fem.functionspace(domain, ("DG", 0))
        w = ufl.TestFunction(W)

        # Project f onto V (data oscillation term): f - f_h
        V = uh.function_space
        fh = project_to_space(f_rhs, V, prefix="proj_f_")

        # Cell residual: f_rhs - f_h (computed in V)
        r = f_rhs - fh

        # Cell size h_T (UFL: CellDiameter)
        h = ufl.CellDiameter(domain)

        # Flux jump term on interior facets. For CG1, the cell residual alone can be
        # too weak because Δu_h is zero inside each cell, while flux jumps still capture
        # where the solution is under-resolved.
        n = ufl.FacetNormal(domain)
        jump_flux = ufl.jump(kappa * ufl.grad(uh), n)
        facet_term = ufl.avg(h) * (jump_flux**2) * ufl.avg(w) * ufl.dS

        # Assemble: ∫_T h^2 r^2 w dx  +  ∫_F h [[∇u·n]]^2 avg(w) dS
        form = fem.form((h**2) * (r**2) * w * ufl.dx + facet_term)

        vec = fem_petsc.assemble_vector(form)
        vec.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)

        eta_sq = vec.array.copy()
        # Ensure non-negative numerical noise doesn't break sqrt
        eta_sq = np.maximum(eta_sq, 0.0)
        eta = np.sqrt(eta_sq)
        return eta
    return (cell_residual_indicator,)


@app.cell
def _(mo):
    mo.md(r"""
    ## Mark & refine

    We mark the top **fraction** of cells (e.g., 20%) with the largest indicators and refine them.
    """)
    return


@app.cell
def _(mesh, np):
    def mark_cells_for_refinement(domain, eta, fraction=0.2):
        """
        Return a boolean array of length num_cells, True for cells to refine.
        Marks the top `fraction` fraction of cells by eta.
        """
        num_cells = domain.topology.index_map(domain.topology.dim).size_local
        assert eta.shape[0] == num_cells

        if fraction <= 0.0:
            return np.zeros(num_cells, dtype=np.int8)
        if fraction >= 1.0:
            return np.ones(num_cells, dtype=np.int8)

        k = max(1, int(np.ceil(fraction * num_cells)))
        # threshold = k-th largest
        thresh = np.partition(eta, -k)[-k]
        markers = (eta >= thresh).astype(np.int8)
        return markers

    def refine_marked(domain, markers):
        """
        Refine marked cells. `markers` must be int8/bool array per local cell.

        Note: In DOLFINx 0.10, `dolfinx.mesh.refine` refines by splitting *edges*,
        so we first convert marked cells -> incident edges.
        """
        tdim = domain.topology.dim
        cells = np.arange(domain.topology.index_map(tdim).size_local, dtype=np.int32)
        marked_cells = cells[markers.astype(bool)]
        if marked_cells.size == 0:
            return domain

        # Build cell -> edge connectivity and collect incident edges.
        domain.topology.create_entities(1)
        domain.topology.create_connectivity(tdim, 1)
        cell_to_edges = domain.topology.connectivity(tdim, 1)
        edges = np.unique(np.hstack([cell_to_edges.links(c) for c in marked_cells])).astype(np.int32)

        refined, _parent_cell, _parent_facet = mesh.refine(domain, edges=edges)
        return refined
    return mark_cells_for_refinement, refine_marked


@app.cell
def _(dplot, np):
    def triangulation_from_dolfinx(obj):
        """
        Build (x, y, triangles) for matplotlib.tri.Triangulation from a 2D dolfinx mesh
        or from a Lagrange function space (via `dolfinx.plot.vtk_mesh`).
        """
        topology, _cell_types, geometry = dplot.vtk_mesh(obj)
        num_nodes = int(topology[0])
        topo = topology.reshape(-1, num_nodes + 1)
        if not np.all(topo[:, 0] == num_nodes):
            raise ValueError("Unexpected VTK topology encoding.")
        cells = topo[:, 1:]
        xy = geometry[:, :2]
        return xy[:, 0], xy[:, 1], cells
    return (triangulation_from_dolfinx,)


@app.cell
def _():
    import matplotlib.tri as mtri
    return (mtri,)


@app.cell
def _():
    import pyvista as pv
    return (pv,)


@app.cell
def _():
    import importlib.util

    trame_vtk_available = importlib.util.find_spec("trame_vtk") is not None
    return (trame_vtk_available,)


@app.cell
def _():
    import trimesh
    from trimesh.viewer.notebook import scene_to_notebook
    from trimesh.visual import color as tm_color
    return trimesh, scene_to_notebook, tm_color


@app.cell
def _(mtri, plt, triangulation_from_dolfinx):
    def plot_mesh(
        domain,
        eta=None,
        title="Mesh",
        cmap="magma",
        label=r"cell indicator $\eta_T$",
        use_log=True,
    ):
        x, y, cells = triangulation_from_dolfinx(domain)
        tri = mtri.Triangulation(x, y, cells)
        fig, ax = plt.subplots()
        ax.triplot(tri, linewidth=0.5, color="tab:blue", alpha=0.65)

        if eta is not None:
            from matplotlib.colors import LogNorm, Normalize

            # `eta` is a per-cell quantity; show it as a light cell fill.
            eta = eta[: cells.shape[0]]
            eta_pos = eta[eta > 0]
            if use_log and eta_pos.size > 0 and float(eta_pos.min()) < float(eta.max()):
                norm = LogNorm(vmin=float(eta_pos.min()), vmax=float(eta.max()))
            else:
                norm = Normalize(vmin=float(eta.min()), vmax=float(eta.max()))

            tpc = ax.tripcolor(
                tri, facecolors=eta, shading="flat", cmap=cmap, norm=norm, alpha=0.45
            )

            fig.colorbar(tpc, ax=ax, label=label)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.tight_layout()
        return fig
    return (plot_mesh,)


@app.cell
def _(mtri, plt, triangulation_from_dolfinx):
    def plot_solution(domain, uh, title="Solution u_h"):
        # Use dof coordinates (not raw mesh vertices) to match `uh.x.array`.
        x, y, cells = triangulation_from_dolfinx(uh.function_space)
        tri = mtri.Triangulation(x, y, cells)

        uvals = uh.x.array.real

        fig, ax = plt.subplots()
        tpc = ax.tripcolor(tri, uvals, shading="gouraud")
        fig.colorbar(tpc, ax=ax)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.tight_layout()
        return fig
    return (plot_solution,)


@app.cell
def _(build_problem, domain0, fem, plot_mesh, project_to_space):
    # Plot the diffusion coefficient on the initial mesh.
    _V0, _kappa_expr0, _kappa0, _f0, _a0, _L0, _bcs0 = build_problem(domain0, degree=1)
    W0 = fem.functionspace(domain0, ("DG", 0))
    kappa_dg0 = project_to_space(_kappa0, W0, prefix="proj_kappa_dg0_")

    fig_f = plot_mesh(
        domain0,
        eta=kappa_dg0.x.array.real,
        title=r"Diffusion coefficient $\kappa(x)$ (DG0 projection on initial mesh)",
        cmap="magma",
        label=r"$\kappa$",
        use_log=True,
    )
    fig_f
    return (fig_f,)


@app.cell
def _(mo):
    mo.md(r"""
    ## Run AMR loop

    We track:
    - number of cells,
    - sum of indicators \(\sum_T \eta_T\) (a rough proxy for "estimated error"),
    - and visualize the mesh concentrating near the source.
    """)
    return


@app.cell
def _(
    MPI,
    cell_residual_indicator,
    comm,
    mark_cells_for_refinement,
    np,
    refine_marked,
    solve_poisson,
):
    def amr_solve(domain, nsteps=4, fraction=0.2):
        history = []
        current = domain
        for k in range(nsteps):
            V, kappa_expr, kappa, f_rhs, uh = solve_poisson(current)
            eta = cell_residual_indicator(current, kappa, f_rhs, uh)
            num_cells_local = current.topology.index_map(current.topology.dim).size_local
            num_cells = comm.allreduce(num_cells_local, op=MPI.SUM)
            eta_sum = comm.allreduce(float(np.sum(eta)), op=MPI.SUM)

            history.append(
                dict(
                    iter=k,
                    num_cells=num_cells,
                    eta_sum=eta_sum,
                    domain=current,
                    uh=uh,
                    eta=eta,
                    kappa=kappa,
                    kappa_expr=kappa_expr,
                )
            )

            markers = mark_cells_for_refinement(current, eta, fraction=fraction)
            current = refine_marked(current, markers)

        return history
    return (amr_solve,)


@app.cell
def _(mo):
    nsteps = mo.ui.slider(1, 6, value=4, label="AMR steps")
    fraction = mo.ui.slider(0.05, 0.5, value=0.2, step=0.05, label="Mark fraction")

    return fraction, nsteps


@app.cell
def _(amr_solve, domain0, fraction, mo, nsteps):
    history = amr_solve(domain0, nsteps=int(nsteps.value), fraction=float(fraction.value))
    controls = mo.hstack([nsteps, fraction])
    controls
    return (history,)


@app.cell
def _(history, mo, plt, rank):
    if rank != 0:
        mo.md("Running under MPI: plots shown only on rank 0.")
        raise SystemExit

    iters = [h["iter"] for h in history]
    ncells = [h["num_cells"] for h in history]
    eta_sum = [h["eta_sum"] for h in history]

    fig, ax = plt.subplots()
    ax.plot(iters, ncells, marker="o")
    ax.set_xlabel("AMR iteration")
    ax.set_ylabel("Number of cells")
    ax.set_title("Mesh growth under AMR")
    fig.tight_layout()

    fig2, ax2 = plt.subplots()
    ax2.plot(iters, eta_sum, marker="o")
    ax2.set_xlabel("AMR iteration")
    ax2.set_ylabel(r"$\sum_T \eta_T$")
    ax2.set_title("Indicator sum (rough error proxy)")
    fig2.tight_layout()
    return fig, fig2


@app.cell
def _(fig, fig2, mo):
    def _render():
        mo.md("### Diagnostics")
        return mo.vstack([fig, fig2])

    _view = _render()
    _view
    return (_view,)


@app.cell
def _(history, mo, rank):
    if rank != 0:
        raise SystemExit

    selector = mo.ui.slider(0, len(history) - 1, value=len(history) - 1, label="Show iteration")
    return (selector,)


@app.cell
def _(history, mo, plot_mesh, plot_solution, rank, selector):
    def _render():
        if rank != 0:
            raise SystemExit

        i = int(selector.value)
        dom = history[i]["domain"]
        uh = history[i]["uh"]
        eta = history[i].get("eta")

        figm = plot_mesh(dom, eta=eta, title=f"Mesh + indicator after iteration {i}")
        figs = plot_solution(dom, uh, title=f"u_h after iteration {i}")

        return mo.vstack([selector, figm, figs])

    _view = _render()
    _view
    return (_view,)


@app.cell
def _(
    dplot,
    history,
    mo,
    np,
    pv,
    rank,
    scene_to_notebook,
    selector,
    tm_color,
    trame_vtk_available,
    trimesh,
):
    def _render():
        if rank != 0:
            raise SystemExit

        i = int(selector.value)
        uh = history[i]["uh"]

        # Build a surface mesh using dof coordinates so point_data aligns with uh.x.array.
        topology, cell_types, geometry = dplot.vtk_mesh(uh.function_space)
        points = np.zeros((geometry.shape[0], 3), dtype=float)
        points[:, :2] = geometry[:, :2]
        points[:, 2] = 3.0 * uh.x.array.real

        # VTK encodes each cell as: [num_nodes, n0, n1, n2]; drop the leading
        # count to get triangle indices, then orient them so normals point +z.
        num_nodes = int(topology[0])
        topo = topology.reshape(-1, num_nodes + 1)
        faces_tri = topo[:, 1:].astype(np.int64)
        v1 = points[faces_tri[:, 1]] - points[faces_tri[:, 0]]
        v2 = points[faces_tri[:, 2]] - points[faces_tri[:, 0]]
        normals = np.cross(v1, v2)
        flip = normals[:, 2] < 0
        faces_oriented = faces_tri.copy()
        faces_oriented[flip] = faces_oriented[flip][:, [0, 2, 1]]

        if trame_vtk_available:
            # Build a PyVista surface and embed it via exported HTML (off-screen).
            pv.OFF_SCREEN = True
            pl = pv.Plotter(off_screen=True)

            # PyVista expects faces in VTK-style flattened format: [3, i0, i1, i2, 3, ...]
            faces_vtk = np.hstack(
                [np.full((faces_oriented.shape[0], 1), 3, dtype=np.int64), faces_oriented]
            ).ravel()

            surface = pv.PolyData(points, faces_vtk)
            surface.point_data["u"] = uh.x.array.real

            pl.add_mesh(surface, scalars="u", cmap="viridis", show_edges=False)
            pl.add_axes()
            pl.view_isometric()
            pl.set_background("white")

            from pathlib import Path
            import base64

            html_path = (
                Path(__file__).resolve().parent / ".cache" / f"pyvista_surface_iter_{i}.html"
            )
            html_path.parent.mkdir(parents=True, exist_ok=True)
            pl.export_html(str(html_path))
            html = html_path.read_text(encoding="utf-8")
            html_b64 = base64.b64encode(html.encode("utf-8")).decode("ascii")
            iframe = (
                f'<iframe src="data:text/html;base64,{html_b64}" '
                f'style="width: 100%; height: 600px; border: 0;"></iframe>'
            )

            class _HtmlRepr:
                def __init__(self, html: str):
                    self._html = html

                def _repr_html_(self) -> str:
                    return self._html

            return mo.vstack(
                [mo.md("### 3D surface view (rotate with mouse, PyVista)"), _HtmlRepr(iframe)]
            )

        # Fallback: embed via trimesh notebook viewer if trame-vtk is missing.
        mesh_tm = trimesh.Trimesh(vertices=points, faces=faces_oriented, process=True)
        mesh_tm.merge_vertices()
        colors = tm_color.interpolate(uh.x.array.real, color_map="viridis")
        mesh_tm.visual.vertex_colors = colors
        scene = trimesh.Scene(mesh_tm)
        view = scene_to_notebook(scene, height=600)
        return mo.vstack([mo.md("### 3D surface view (rotate with mouse, trimesh)"), view])

    _pv_view = _render()
    _pv_view
    return (_pv_view,)


@app.cell
def _(mo):
    mo.md(r"""
    ## What to emphasize in class

    - Even a simple indicator captures where the solution is difficult (near the source).
    - AMR reallocates DOFs **where they matter**, usually outperforming uniform refinement.
    - The full story includes:
      - flux-jump terms on facets,
      - goal-oriented (DWR) indicators,
      - coarsening and refinement strategies,
      - and adaptivity for time-dependent PDEs.
    """)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
# Extension to MKP: Nitsche Method via an Extended-Domain Interpretation

**Goal**

Derive Nitsche’s method *from the honest split-domain calculation*:

1) Start from a **two-domain** (interior + exterior) weak statement with the **full flux jump** on the interface.  
2) Impose the model “**outside equals prescribed data**” and eliminate the exterior **by calculation** (no reverse engineering).  
3) Identify what remains missing (a closure on the interface) and introduce the **local, symmetric weak enforcement** that yields Nitsche.

This makes the contrast to MKP clear:
- MKP-type approaches enforce boundary conditions by adding a **forcing/penalty source** (often volumetric or delta-supported).
- Nitsche enforces boundary conditions by modifying the **variational operator and load** with **interface terms**, consistent with eliminating an exterior domain.

We use a scalar diffusion model for clarity; the structure generalizes to Stokes/Navier–Stokes.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 1. Model problem and extension setup

Physical domain: \(\Omega \subset \mathbb{R}^d\), boundary \(\Gamma\).

Diffusion model:

\[
-\nabla\cdot(\kappa \nabla u)=f \quad \text{in }\Omega,
\qquad
u=g \quad \text{on }\Gamma.
\]

**Fictitious-domain / immersed viewpoint**

Introduce an exterior region \(\Omega^+\) sharing the interface \(\Gamma\) with the interior:

\[
\Omega^- := \Omega,\qquad \Gamma = \partial\Omega^- \cap \partial\Omega^+.
\]

Choose the unit normal \(n\) pointing **from \(\Omega^-\) to \(\Omega^+\)**.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 2. Strong equations in each subdomain

Interior equation:

\[
-\nabla\cdot(\kappa^- \nabla u^-)=f \quad \text{in }\Omega^-.
\]

Exterior modeling assumption (prescribed extension):

\[
u^+ = g \quad \text{in }\Omega^+.
\]

If the diffusion operator is taken to hold in \(\Omega^+\) as well, then prescribing \(u^+=g\) implies an exterior forcing:

\[
-\nabla\cdot(\kappa^+\nabla u^+) = -\nabla\cdot(\kappa^+\nabla g) =: f^+ \quad \text{in }\Omega^+.
\]

This “implied \(f^+\)” is the clean way to see how a prescribed exterior state can correspond to a *classical source*—but only in the exterior region.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 3. Weak statement on the combined domain

Let \(v\) be a test function defined on \(\Omega^- \cup \Omega^+\) (not restricted to vanish on \(\Gamma\)).

Multiply by \(v\) and integrate:

\[
\int_{\Omega^-} \big(-\nabla\cdot(\kappa^- \nabla u^-)\big) v
+
\int_{\Omega^+} \big(-\nabla\cdot(\kappa^+ \nabla u^+)\big) v
=
\int_{\Omega^-} f v + \int_{\Omega^+} f^+ v.
\]

Integrate by parts **separately** in each subdomain and then add the two identities.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 4. Integration by parts and the full flux jump

After integration by parts in \(\Omega^-\) and \(\Omega^+\) and adding, the interface term combines into a **jump of normal flux**:

\[
\int_{\Omega^-}\kappa^- \nabla u^- \cdot \nabla v
+
\int_{\Omega^+}\kappa^+ \nabla u^+ \cdot \nabla v
-
\int_{\Gamma} [\![ \kappa \nabla u \cdot n ]\!]\, v
=
\int_{\Omega^-} f v + \int_{\Omega^+} f^+ v.
\]

With the chosen normal \(n\) (out of \(\Omega^-\)):

\[
[\![ \kappa \nabla u \cdot n ]\!] :=
(\kappa^- \nabla u^- \cdot n) - (\kappa^+ \nabla u^+ \cdot n).
\]

This is the *complete* interface term: the flux from both sides appears automatically from the calculus.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 5. Substitute the prescribed exterior state

Using \(u^+=g\) in \(\Omega^+\):

\[
\nabla u^+ = \nabla g,\qquad (\kappa \nabla u\cdot n)^+ = \kappa^+ \nabla g\cdot n.
\]

The combined-domain identity becomes:

- Interior energy term: \(\int_{\Omega^-} \kappa^- \nabla u^- \cdot \nabla v\)
- Exterior energy term: \(\int_{\Omega^+} \kappa^+ \nabla g \cdot \nabla v\)
- Interface term: \(-\int_\Gamma \big((\kappa^- \nabla u^- \cdot n) - (\kappa^+ \nabla g \cdot n)\big) v\)
- Right-hand side: \(\int_{\Omega^-} f v + \int_{\Omega^+} f^+ v\)

At this point the only “unknown” is \(u^-\) in \(\Omega^-\); everything in \(\Omega^+\) is data.
The next step is to eliminate the exterior **volume** term by calculation.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 6. Eliminate the exterior volume term (the “opposite per partes” step)

Consider the exterior volume term:

\[
\int_{\Omega^+} \kappa^+ \nabla g \cdot \nabla v.
\]

Integrate by parts in \(\Omega^+\):

\[
\int_{\Omega^+} \kappa^+ \nabla g \cdot \nabla v
=
-\int_{\Omega^+} \nabla\cdot(\kappa^+ \nabla g)\, v
+ \int_{\partial\Omega^+} (\kappa^+\nabla g \cdot n^+)\, v.
\]

Using \(f^+ = -\nabla\cdot(\kappa^+\nabla g)\), the exterior **volume** contribution cancels exactly with \(\int_{\Omega^+} f^+ v\).

What remains from the exterior is **boundary/interface data** (plus possibly the outer boundary of \(\Omega^+\); here we assume either it is absent or handled separately).
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 7. After eliminating \(\Omega^+\): interior-only identity

After cancellation, the two-domain formulation reduces to an interior integration-by-parts identity on \(\Omega=\Omega^-\):

\[
\int_{\Omega} \kappa \nabla u \cdot \nabla v
-
\int_{\Gamma} \kappa(\nabla u\cdot n)\, v
=
\int_{\Omega} f v.
\]

This shows what elimination *actually* does:
- The exterior volume physics collapses to interface information.
- Dropping \(\Omega^+\) is legitimate only if we still encode the condition “outside equals \(g\)” through an interface closure.

The remaining task is: impose \(u=g\) on \(\Gamma\) **weakly and stably** without constraining the trial/test space.
That closure is Nitsche’s method.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 8. Interface closure and why Nitsche has the terms it has

Exact elimination of \(\Omega^+\) produces a (generally nonlocal) relation between the interface trace and flux:
a Dirichlet-to-Neumann type operator. This is not practical for most discretizations.

Nitsche replaces that exact exterior response by a **local interface law** that:
- enforces the trace condition \(u=g\) weakly,
- keeps the formulation consistent,
- provides stability/coercivity via a scaled stabilization.

A convenient local closure is Robin-like:

\[
\kappa \nabla u \cdot n + \frac{\gamma \kappa}{h}(u-g) \approx \text{(exterior flux data)}.
\]

To enforce this without losing symmetry and accuracy, the weak form uses:
- flux–trace pairing,
- the symmetric “transpose” pairing (adjoint consistency),
- and stabilization scaled as \(\kappa/h\).

This leads to the standard symmetric Nitsche formulation.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 9. Nitsche formulation: bilinear and linear forms

Write the method as: find \(u_h\) such that \(a(u_h,v_h)=L(v_h)\) for all \(v_h\).

### Bilinear form

\[
a(u,v)
=
\int_{\Omega} \kappa \nabla u\cdot \nabla v \, d\Omega
-\int_{\Gamma} \kappa(\nabla u\cdot n)\, v \, d\Gamma
-\int_{\Gamma} \kappa(\nabla v\cdot n)\, u \, d\Gamma
+\int_{\Gamma} \frac{\gamma \kappa}{h}\, u v \, d\Gamma.
\]

### Linear form

\[
L(v)
=
\int_{\Omega} f v \, d\Omega
-\int_{\Gamma} \kappa(\nabla v\cdot n)\, g \, d\Gamma
+\int_{\Gamma} \frac{\gamma \kappa}{h}\, g v \, d\Gamma.
\]

Interpretation in the extended-domain picture:
- the interior PDE remains unchanged,
- boundary enforcement is encoded through symmetric interface terms,
- stabilization controls coercivity without an artificial volumetric source.
"""
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
## 10. MKP contrast (kept to one clean statement)

In MKP-style enforcement, the boundary condition is often imposed by adding a forcing/penalty term to the equation (or its weak form) that behaves like a **source** driving \(u\to g\).

In the extended-domain interpretation above, Nitsche instead corresponds to:
- eliminating the exterior region,
- then approximating the (otherwise nonlocal) interface response by **boundary-only variational couplings**.

So the “extra effect” shows up as *interface operator + interface load*, not as a volumetric source in \(\Omega\).
"""
    )
    return


if __name__ == "__main__":
    app.run()
