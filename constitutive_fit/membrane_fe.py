"""
Small finite element solver for the incompressible plane-stress membrane (constitutive.py) on a
structured grid of bilinear (Q1) elements, with prescribed displacements on the whole boundary.

Used to generate synthetic data for the Stage 2 virtual-field equations (stage2_synthetic.py): the
nodes are the DIC grid points, so a solution can be passed through the same F pipeline as the
measured tracks.

    mesh = Mesh(Xg)                        # Xg (ny, nx, 2) reference node coordinates (mm), regular grid
    x, info = mesh.solve(xb, alpha_gp, P, H, x0)
        xb: (ny, nx, 2) positions; only the boundary nodes are used (prescribed)
        alpha_gp: (ne_r, ne_c, 4) stiffness factor at each Gauss point (scales Cgr, C1_el, C1_smc, C1_col), or
                  (ne_r, ne_c, 4, n) with P a list of n parameter vectors whose energies are summed
        x0: initial guess for all nodes (boundary values are overwritten by xb)
    mesh.reactions(x, alpha_gp, P, H)  -> (ny, nx, 2) nodal forces (mN) the boundary must apply

Energy: Pi = sum_e sum_g H alpha_g W(F_g) A_e / 4 (2 x 2 Gauss quadrature; W from constitutive.strain_energy,
in kPa = mN/mm^2; H in mm), minimized over the free nodes with L-BFGS and the analytic gradient
dPi/dx_a = sum_g alpha_g N(F_g) . grad N_a(g) A_e / 4, N = first_pk_membrane(F, P, H).
"""
import numpy as np
from scipy.optimize import minimize

import constitutive as M

_G = np.array([-1, 1]) / np.sqrt(3)                       # 2-point Gauss abscissae on [-1, 1]


class Mesh:
    def __init__(self, Xg):
        self.Xg = np.asarray(Xg, float)
        ny, nx = self.Xg.shape[:2]
        self.ny, self.nx = ny, nx
        hx = np.diff(self.Xg[0, :, 0])
        hy = np.diff(self.Xg[:, 0, 1])
        if not (np.allclose(hx, hx[0]) and np.allclose(hy, hy[0])):
            raise ValueError("Mesh expects a regular grid")
        self.hx, self.hy = hx[0], hy[0]                   # element size (hy may be negative; kept signed)
        self.area = abs(self.hx * self.hy)
        # shape-function gradients w.r.t. reference coords at the 4 Gauss points, corners (bl, br, tr, tl)
        # in grid-index order: corner 0 = (r, c), 1 = (r, c+1), 2 = (r+1, c+1), 3 = (r+1, c)
        dN = []
        for gy in _G:
            for gx in _G:
                s, t = (gx + 1) / 2, (gy + 1) / 2                 # local coords in [0, 1]
                dNds = np.array([-(1 - t), (1 - t), t, -t])
                dNdt = np.array([-(1 - s), -s, s, (1 - s)])
                dN.append(np.stack([dNds / self.hx, dNdt / self.hy], 1))   # (4 corners, 2)
        self.dN = np.array(dN)                                    # (4 gp, 4 corners, 2)
        free = np.ones((ny, nx), bool)
        free[0], free[-1], free[:, 0], free[:, -1] = False, False, False, False
        self.free = free
        # Gauss point reference coordinates (ne_r, ne_c, 4, 2), for assigning material
        c = self._corners(self.Xg)
        gp = []
        for gy in _G:
            for gx in _G:
                s, t = (gx + 1) / 2, (gy + 1) / 2
                Nn = np.array([(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t])
                gp.append(np.einsum("a,rcad->rcd", Nn, c))
        self.gauss_points = np.stack(gp, 2)

    @staticmethod
    def _corners(x):
        return np.stack([x[:-1, :-1], x[:-1, 1:], x[1:, 1:], x[1:, :-1]], 2)   # (ne_r, ne_c, 4, 2)

    def F(self, x):
        """Deformation gradient at every Gauss point: (ne_r, ne_c, 4, 2, 2)."""
        return np.einsum("rcai,gaJ->rcgiJ", self._corners(x), self.dN)

    @staticmethod
    def _parts(alpha_gp, P):
        """Material as a sum of parts: a single (alpha_gp, P), or alpha_gp (..., n) with P a list of n parameter
        vectors (e.g. matrix/elastin and fiber constants scaled separately, ladder level (b))."""
        if isinstance(P, (list, tuple)):
            return [(alpha_gp[..., c], np.asarray(Pc, float)) for c, Pc in enumerate(P)]
        return [(alpha_gp, P)]

    def energy_grad(self, x, alpha_gp, P, H):
        F = self.F(x)
        if np.any(np.linalg.det(F) <= 0):               # inverted element: large finite value keeps line searches sane
            return 1e30, np.zeros_like(x)
        Pi, N = 0.0, 0.0
        for al, Pc in self._parts(alpha_gp, P):
            w = H * al * self.area / 4
            Pi = Pi + np.sum(w * M.strain_energy(F, Pc))
            N = N + M.first_pk_membrane(F, Pc, 1.0) * w[..., None, None]
        fe = np.einsum("rcgiJ,gaJ->rcai", N, self.dN)            # element nodal forces (ne_r, ne_c, 4, 2)
        g = np.zeros_like(x)
        g[:-1, :-1] += fe[:, :, 0]; g[:-1, 1:] += fe[:, :, 1]
        g[1:, 1:] += fe[:, :, 2]; g[1:, :-1] += fe[:, :, 3]
        return Pi, g

    def solve(self, xb, alpha_gp, P, H, x0, gtol=1e-6, maxiter=20000, scale=0.01):
        """Minimize the energy over the free nodes. The unknowns are scaled (x = x0 + scale * v, mm) so that
        L-BFGS-B's first trial step, of unit length, moves nodes by about `scale` mm rather than 1 mm."""
        x = np.array(x0, float)
        bnd = ~self.free
        x[bnd] = xb[bnd]
        xf0 = x[self.free].copy()

        def fun(v):
            x[self.free] = xf0 + scale * v.reshape(-1, 2)
            Pi, g = self.energy_grad(x, alpha_gp, P, H)
            return Pi, scale * g[self.free].ravel()
        r = minimize(fun, np.zeros(xf0.size), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=maxiter, maxcor=30, gtol=gtol * scale, ftol=1e-15))
        x[self.free] = xf0 + scale * r.x.reshape(-1, 2)
        _, g = self.energy_grad(x, alpha_gp, P, H)
        reac = g[bnd]
        info = dict(success=bool(r.success), nit=int(r.nit), message=str(r.message),
                    max_free_force=float(np.abs(g[self.free]).max()),
                    typical_boundary_force=float(np.median(np.abs(reac[np.abs(reac) > 0])) if np.any(reac) else 0))
        return x, info

    def reactions(self, x, alpha_gp, P, H):
        """Nodal forces (mN) the supports must apply; nonzero only on the boundary at equilibrium."""
        _, g = self.energy_grad(x, alpha_gp, P, H)
        return g
