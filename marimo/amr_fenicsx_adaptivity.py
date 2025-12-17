# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "marimo==0.18.4",
#   "numpy",
#   "matplotlib",
#   "mpi4py",
#   "petsc4py",
#   "fenics-dolfinx==0.10.0",
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
    import numpy as np
    import matplotlib.pyplot as plt

    from mpi4py import MPI
    from petsc4py import PETSc

    import ufl
    from dolfinx import fem, mesh
    return MPI, PETSc, fem, mesh, np, plt, ufl


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
    -\Delta u = f \quad \text{in } \Omega=(0,1)^2, \qquad
    u=0 \quad \text{on } \partial \Omega.
    \]

    with a Gaussian right-hand side centered at \((0.7, 0.35)\).
    """)
    return


@app.cell
def _(PETSc, fem, mesh, np, ufl):
    def build_problem(domain: mesh.Mesh, degree: int = 1):
        """Return (V, u, f, a, L, bcs) for Poisson with homogeneous Dirichlet BC."""
        V = fem.functionspace(domain, ("Lagrange", degree))

        u = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        x = ufl.SpatialCoordinate(domain)
        # A localized source: makes a sharp feature that AMR will chase.
        x0, y0, alpha = 0.7, 0.35, 200.0
        f_expr = ufl.exp(-alpha * ((x[0] - x0) ** 2 + (x[1] - y0) ** 2))

        a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
        L = f_expr * v * ufl.dx

        # Homogeneous Dirichlet BC on the whole boundary.
        def boundary(x):
            # "near" is not required; for unit square with float coords, this is ok.
            return np.isclose(x[0], 0.0) | np.isclose(x[0], 1.0) | np.isclose(x[1], 0.0) | np.isclose(x[1], 1.0)

        facets = mesh.locate_entities_boundary(domain, domain.topology.dim - 1, boundary)
        dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, facets)
        bc = fem.dirichletbc(PETSc.ScalarType(0), dofs, V)

        return V, f_expr, a, L, [bc]
    return (build_problem,)


@app.cell
def _(comm, mesh):
    # Start with a coarse mesh.
    domain0 = mesh.create_unit_square(comm, 10, 10, cell_type=mesh.CellType.triangle)
    return (domain0,)


@app.cell
def _(build_problem, fem):
    def solve_poisson(domain):
        V, f_expr, a, L, bcs = build_problem(domain, degree=1)

        problem = fem.petsc.LinearProblem(
            a, L, bcs=bcs,
            petsc_options={
                "ksp_type": "cg",
                "pc_type": "hypre",
                "ksp_rtol": 1.0e-10,
            },
        )
        uh = problem.solve()
        uh.name = "u"
        return V, f_expr, uh
    return (solve_poisson,)


@app.cell
def _(mo):
    mo.md(r"""
    ## Residual-based local indicator (teaching version)

    For each cell \(T\), define

    \[
    \eta_T^2 \approx h_T^2 \| f + \Delta u_h \|^2_{L^2(T)}
    \]

    where \(h_T\) is a cell size measure.

    This is a simplified residual indicator. In a full estimator you would also
    include **edge jump terms** of the flux to improve reliability, but this is
    enough to demonstrate adaptive refinement behavior in 90 minutes.
    """)
    return


@app.cell
def _(PETSc, fem, np, ufl):
    def cell_residual_indicator(domain, f_expr, uh):
        """
        Compute a per-cell indicator eta_T using a simple cell residual:
            r = f + div(grad(uh))   (note: for Poisson with -Δu=f => Δu = -f)
        We assemble ∫_T h_T^2 r^2 dx for each cell T by using a DG0 test space.
        """
        # DG0 space: one DOF per cell
        W = fem.functionspace(domain, ("DG", 0))
        w = ufl.TestFunction(W)

        # Residual for -Δu = f is: r = f + Δu_h
        r = f_expr + ufl.div(ufl.grad(uh))

        # Cell size h_T (UFL: CellDiameter)
        h = ufl.CellDiameter(domain)

        # We assemble the vector with entries: ∫ h^2 r^2 * w dx
        # Since w is DG0 test function, this produces one value per cell.
        form = fem.form((h**2) * (r**2) * w * ufl.dx)

        vec = fem.petsc.assemble_vector(form)
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
        """
        # dolfinx.mesh.refine accepts a meshtag or boolean array in some versions.
        # Use a MeshTags for portability.
        tdim = domain.topology.dim
        cells = np.arange(domain.topology.index_map(tdim).size_local, dtype=np.int32)
        marked_cells = cells[markers.astype(bool)]
        # Create meshtags with marked cells = 1
        mt = mesh.meshtags(domain, tdim, marked_cells, np.ones_like(marked_cells, dtype=np.int32))
        refined = mesh.refine(domain, mt)
        return refined
    return mark_cells_for_refinement, refine_marked


@app.cell
def _(mesh):
    def triangulation_from_dolfinx(domain: mesh.Mesh):
        """
        Build (x, y, triangles) for matplotlib.tri.Triangulation from a 2D dolfinx mesh.
        Works for triangular meshes.
        """
        tdim = domain.topology.dim
        domain.topology.create_connectivity(tdim, 0)
        conn = domain.topology.connectivity(tdim, 0)
        cells = conn.array.reshape(-1, conn.num_links(0))  # triangles
        xy = domain.geometry.x[:, :2]
        return xy[:, 0], xy[:, 1], cells
    return (triangulation_from_dolfinx,)


@app.cell
def _(plt, triangulation_from_dolfinx):
    import matplotlib.tri as mtri

    def plot_mesh(domain, title="Mesh"):
        x, y, cells = triangulation_from_dolfinx(domain)
        tri = mtri.Triangulation(x, y, cells)
        fig, ax = plt.subplots()
        ax.triplot(tri, linewidth=0.4)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.tight_layout()
        return fig
    return (plot_mesh,)


@app.cell
def _(plt, triangulation_from_dolfinx):
    import matplotlib.tri as mtri

    def plot_solution(domain, uh, title="Solution u_h"):
        x, y, cells = triangulation_from_dolfinx(domain)
        tri = mtri.Triangulation(x, y, cells)

        # Evaluate at mesh vertices by pulling the dof values for CG1.
        # For CG1, dofs coincide with vertices on this simple mesh.
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
            V, f_expr, uh = solve_poisson(current)
            eta = cell_residual_indicator(current, f_expr, uh)
            num_cells_local = current.topology.index_map(current.topology.dim).size_local
            num_cells = comm.allreduce(num_cells_local, op=MPI.SUM)
            eta_sum = comm.allreduce(float(np.sum(eta)), op=MPI.SUM)

            history.append(
                dict(iter=k, num_cells=num_cells, eta_sum=eta_sum, domain=current, uh=uh, eta=eta)
            )

            markers = mark_cells_for_refinement(current, eta, fraction=fraction)
            current = refine_marked(current, markers)

        return history
    return (amr_solve,)


@app.cell
def _(amr_solve, domain0, mo):
    nsteps = mo.ui.slider(1, 6, value=4, label="AMR steps")
    fraction = mo.ui.slider(0.05, 0.5, value=0.2, step=0.05, label="Mark fraction")
    history = amr_solve(domain0, nsteps=int(nsteps.value), fraction=float(fraction.value))
    mo.vstack([mo.hstack([nsteps, fraction])])
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
    mo.md("### Diagnostics")
    mo.vstack([mo.pyplot(fig), mo.pyplot(fig2)])
    return


@app.cell
def _(history, mo, plot_mesh, plot_solution, rank):
    if rank != 0:
        raise SystemExit

    selector = mo.ui.slider(0, len(history) - 1, value=len(history) - 1, label="Show iteration")
    i = int(selector.value)
    dom = history[i]["domain"]
    uh = history[i]["uh"]

    figm = plot_mesh(dom, title=f"Mesh after iteration {i}")
    figs = plot_solution(dom, uh, title=f"u_h after iteration {i}")

    mo.vstack([selector, mo.pyplot(figm), mo.pyplot(figs)])
    return


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


if __name__ == "__main__":
    app.run()
