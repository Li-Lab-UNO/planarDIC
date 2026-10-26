"""
Finite element solver for the four-fiber model (constitutive.py) as a 3D solid: a slab on the DIC grid, hexahedral
(Q1) elements, several element layers through the thickness, prescribed displacements on the side faces.

Geometry: in-plane node positions (z, theta) are the DIC grid points (mm, specimen frame); the third coordinate r is
through the thickness. With faces="bottom" the top face (the imaged face) is flat at r = 0 and the bottom face at
r = -H(z, theta); with faces="both" the mid-plane is flat at r = 0 and the faces at +-H/2, so a thickness change sits
on both faces equally and puts no load off the mid-plane. Element layers are equally spaced in each column.

Material: the same strain energy as the membrane model, written for a 3D deformation gradient, with near
incompressibility enforced by a penalty on each element's mean volume change (mean-dilatation, Q1/P0), which avoids
volumetric locking of fully integrated hexahedra:
    Pi = sum_e alpha_e [ sum_g W(F_g) w_g det J0_g  +  V0_e (-Cgr ln Theta_e) ]  +  sum_e V0_e kappa/2 (ln Theta_e)^2
    W(F) = Cgr/2 (I1 - 3) + sum_i C1_i/(4 C2_i) (exp(C2_i <I4_i - 1>^2) - 1)
with Theta_e = current / reference element volume; fibers lie in the z-theta plane (elastin along z, SMC along theta,
collagen at +-gamma from z). At Theta = 1 the -Cgr ln Theta term makes the reference stress free, and for a
homogeneous incompressible deformation the stress equals the membrane model's (checked below). alpha_e scales Cgr and
every C1 (a stiffness factor per element); kappa is a numerical penalty and is not scaled.

Boundary conditions: in-plane positions prescribed on every node of the side faces (all layers: the tines pierce the
full thickness), r prescribed on the side nodes of the middle node layer (removes the out-of-plane rigid motions);
top and bottom faces free.

Solution: Newton's method on the energy with a sparse direct solve; the element tangent is the one-sided finite
difference of the analytic element forces, and a step is converged when the largest free nodal force is below 1e-6
of the largest reaction. From a start that already meets the boundary (a membrane solution, from_membrane), the
through-thickness positions are relaxed first and then all free positions. Otherwise, or if that fails, the boundary
is applied in load steps from the reference, each started with the linearized predictor of the prescribed increment
(backtracked if it inverts elements); steps are halved when Newton fails.

    slab = Slab(Xg, H, layers=4)                    # Xg (ny, nx, 2) in-plane reference grid (mm); H scalar or (ny, nx)
    x, info = slab.solve(xb_inplane, P, alpha=None)  # xb (ny, nx, 2) in-plane positions of the side nodes (mm)
    slab.top(x)                                      # (ny, nx, 3) top-face node positions
    x0 = slab.from_membrane(xm)                      # initial guess from a membrane solution (ny, nx, 2)
    x0 = slab.from_coarser(x2, 2)                    # ... or from a solution with fewer element layers

    .venv/bin/python constitutive_fit/solid_fe.py    runs the self-checks
"""
import time

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

KAPPA_FACTOR = 2000.0          # kappa = this x (Cgr + sum C1): J - 1 ~ stress / kappa
_G = np.array([-1.0, 1.0]) / np.sqrt(3.0)
_LOCAL = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                   [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)   # (xi along z, eta along theta, zeta along r)


def _dN():
    """Shape-function derivatives w.r.t. (xi, eta, zeta) at the 8 Gauss points: (8 gp, 8 nodes, 3)."""
    out = []
    for z in _G:
        for e in _G:
            for x in _G:
                g = np.array([x, e, z])
                d = np.empty((8, 3))
                for a, (xa, ea, za) in enumerate(_LOCAL):
                    d[a] = 0.125 * np.array([xa * (1 + ea * g[1]) * (1 + za * g[2]),
                                             ea * (1 + xa * g[0]) * (1 + za * g[2]),
                                             za * (1 + xa * g[0]) * (1 + ea * g[1])])
                out.append(d)
    return np.array(out)


def fiber_dirs(gamma_deg):
    g = np.radians(gamma_deg)
    return np.array([[1, 0, 0], [0, 1, 0], [np.cos(g), np.sin(g), 0], [np.cos(g), -np.sin(g), 0]], float)


def energy_density(F, P):
    """W (kPa) of the isochoric-free part (no volumetric term) for F (..., 3, 3)."""
    Cgr, C1e, C2e, C1s, C2s, C1c, C2c, g = P
    C = np.einsum("...ki,...kj->...ij", F, F)
    W = 0.5 * Cgr * (np.trace(C, axis1=-2, axis2=-1) - 3)
    for M, c1, c2 in zip(fiber_dirs(g), (C1e, C1s, C1c, C1c), (C2e, C2s, C2c, C2c)):
        x = np.maximum(np.einsum("i,...ij,j->...", M, C, M) - 1, 0)
        W = W + c1 / (4 * c2) * (np.exp(c2 * x ** 2) - 1)
    return W


def first_pk(F, P):
    """dW/dF (kPa) of energy_density: Cgr F + sum_i C1_i <I4-1> exp(C2_i (I4-1)^2) (F M_i) (x) M_i."""
    Cgr, C1e, C2e, C1s, C2s, C1c, C2c, g = P
    out = Cgr * F
    for M, c1, c2 in zip(fiber_dirs(g), (C1e, C1s, C1c, C1c), (C2e, C2s, C2c, C2c)):
        FM = np.einsum("...ij,j->...i", F, M)
        e = np.einsum("...i,...i->...", FM, FM) - 1
        w = c1 * np.maximum(e, 0) * np.exp(c2 * e ** 2)
        out = out + w[..., None, None] * FM[..., :, None] * M[None, :]
    return out


class Slab:
    def __init__(self, Xg, H, layers=4, faces="bottom"):
        Xg = np.asarray(Xg, float)
        self.ny, self.nx = Xg.shape[:2]
        self.L = int(layers)
        if self.L % 2:
            raise ValueError("layers must be even (the middle node layer holds the out-of-plane constraint)")
        Hn = np.broadcast_to(np.asarray(H, float), (self.ny, self.nx))
        self.H = np.array(Hn)
        if faces not in ("bottom", "both"):
            raise ValueError("faces: 'bottom' (top face flat) or 'both' (mid-plane flat)")
        self.faces = faces
        s = np.linspace(-1.0, 0.0, self.L + 1) + (0.5 if faces == "both" else 0.0)   # bottom to top
        X = np.empty((self.L + 1, self.ny, self.nx, 3))
        X[..., 0], X[..., 1] = Xg[None, ..., 0], Xg[None, ..., 1]
        X[..., 2] = s[:, None, None] * Hn[None]
        self.X = X
        self.nn = X[..., 0].size
        idx = np.arange(self.nn).reshape(self.L + 1, self.ny, self.nx)
        k, i, j = np.meshgrid(np.arange(self.L), np.arange(self.ny - 1), np.arange(self.nx - 1), indexing="ij")
        corners = [(0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0), (1, 0, 0), (1, 0, 1), (1, 1, 1), (1, 1, 0)]
        self.conn = np.stack([idx[k + dk, i + di, j + dj].ravel() for dk, di, dj in corners], 1)   # (ne, 8)
        self.ne = len(self.conn)
        self.elem_shape = (self.L, self.ny - 1, self.nx - 1)
        dN = _dN()
        Xe = X.reshape(-1, 3)[self.conn]                            # (ne, 8, 3)
        J0 = np.einsum("eai,gaj->egij", Xe, dN)                      # dX/dxi
        self.detJ0 = np.linalg.det(J0)
        if np.any(self.detJ0 <= 0):
            raise ValueError("inverted reference element")
        self.gradN = np.einsum("gaj,egji->egai", dN, np.linalg.inv(J0))   # dN/dX (ne, 8 gp, 8 nodes, 3)
        self.V0 = self.detJ0.sum(1)                                  # Gauss weights are 1
        # degrees of freedom: side nodes in-plane prescribed; middle-layer side nodes r prescribed
        side = np.zeros((self.ny, self.nx), bool)
        side[0], side[-1], side[:, 0], side[:, -1] = True, True, True, True
        fixed = np.zeros((self.L + 1, self.ny, self.nx, 3), bool)
        fixed[:, side, 0] = True
        fixed[:, side, 1] = True
        fixed[self.L // 2, side, 2] = True
        self.fixed = fixed.ravel()
        self.free = ~self.fixed
        edof = (3 * self.conn[:, :, None] + np.arange(3)).reshape(self.ne, 24)
        self.edof = edof
        self.rows = np.repeat(edof, 24, 1).ravel()
        self.cols = np.tile(edof, (1, 24)).ravel()

    # ------------------------------------------------------------------ element kernels
    def _kappa(self, P):
        return KAPPA_FACTOR * (P[0] + P[1] + P[3] + P[5])

    def element_forces(self, xe, P, alpha, ge=slice(None)):
        """Element nodal forces (n, 8, 3) (mN for kPa and mm) and element energies for element positions xe (n, 8, 3)."""
        gN = self.gradN[ge]
        dJ0 = self.detJ0[ge]
        F = np.einsum("eai,egaj->egij", xe, gN)                      # (n, 8 gp, 3, 3)
        J = np.linalg.det(F)
        if np.any(J <= 0):
            return None, None
        Finv_T = np.swapaxes(np.linalg.inv(F), -1, -2)
        V0 = self.V0[ge]
        theta = (J * dJ0).sum(1) / V0
        lt = np.log(theta)
        kap = self._kappa(P)
        Cgr = P[0]
        a = alpha[ge]
        Pk = first_pk(F, P) * (a[:, None, None, None] * dJ0[..., None, None])
        dU = (a * (-Cgr) + kap * lt) / theta                         # dU/dTheta per unit reference volume
        Pv = (dU[:, None, None, None] * J[..., None, None] * Finv_T * dJ0[..., None, None])  # d(V0 U)/dF per gp
        f = np.einsum("egij,egaj->eai", Pk + Pv, gN)
        Wd = (energy_density(F, P) * dJ0).sum(1) * a + V0 * (a * (-Cgr) * lt + 0.5 * kap * lt ** 2)
        return f, Wd

    def energy_forces(self, x, P, alpha):
        xe = x[self.conn]
        f, Wd = self.element_forces(xe, P, alpha)
        if f is None:
            return np.inf, None
        g = np.zeros((self.nn, 3))
        np.add.at(g, self.conn, f)
        return float(Wd.sum()), g.ravel()

    def tangent(self, x, P, alpha, h=1e-7, f0=None):
        """Element stiffness by one-sided finite differences of the analytic element forces, symmetrized, assembled."""
        xe = x[self.conn]
        if f0 is None:
            f0, _ = self.element_forces(xe, P, alpha)
        Ke = np.empty((self.ne, 24, 24))
        for k in range(24):
            a, i = divmod(k, 3)
            xp = xe.copy()
            xp[:, a, i] += h
            fp, _ = self.element_forces(xp, P, alpha)
            Ke[:, :, k] = (fp - f0).reshape(self.ne, 24) / h
        Ke = 0.5 * (Ke + np.swapaxes(Ke, 1, 2))
        return sp.coo_matrix((Ke.ravel(), (self.rows, self.cols)), shape=(3 * self.nn, 3 * self.nn)).tocsr()

    # ------------------------------------------------------------------ solution
    def prescribed(self, xb_inplane, x=None):
        """Node positions with the side nodes' in-plane coordinates from xb_inplane (ny, nx, 2) (r unchanged)."""
        out = (self.X if x is None else x.reshape(self.X.shape)).copy()
        side = self.fixed.reshape(self.X.shape)[..., 0]
        out[..., :2] = np.where(side[..., None], np.asarray(xb_inplane)[None], out[..., :2])
        return out.reshape(-1, 3)

    def solve(self, xb_inplane, P, alpha=None, x0=None, steps=4, tol=1e-6, max_iter=20, verbose=False):
        """Equilibrium with the side nodes at xb_inplane, from x0 (default: the reference). Returns (x (nn, 3), info)."""
        P = np.asarray(P, float)
        alpha = np.ones(self.ne) if alpha is None else np.asarray(alpha, float).ravel()
        x = (self.X.reshape(-1, 3) if x0 is None else np.asarray(x0, float).reshape(-1, 3)).copy()
        target = self.prescribed(xb_inplane, x).ravel()
        start = x.ravel().copy()
        fr, fx = self.free, self.fixed
        t0, n_newton, n_steps = time.time(), 0, 0
        if x0 is not None and np.allclose(target, start):
            # the start already meets the boundary (from_membrane): through-thickness positions first, then all
            u = start.copy()
            rdof = np.zeros((self.nn, 3), bool)
            rdof[:, 2] = True
            ok_r, it_r = self._newton(u, P, alpha, tol, max_iter, dofs=rdof.ravel() & fr)
            ok, it = self._newton(u, P, alpha, tol, 2 * max_iter)
            if verbose:
                print(f"    from the start: {it_r} through-thickness + {it} Newton iterations ({time.time() - t0:.1f} s)",
                      flush=True)
            if ok:
                _, R = self.energy_forces(u.reshape(-1, 3), P, alpha)
                return u.reshape(-1, 3), dict(success=True, load=1.0, newton=it_r + it, steps=0, from_start=True,
                                              seconds=time.time() - t0, max_free_force=float(np.abs(R[fr]).max()),
                                              typical_reaction=float(np.median(np.abs(R[fx][np.abs(R[fx]) > 0]))))
            x = self.X.reshape(-1, 3).copy()                         # fall back to load steps from the reference
            start = x.ravel().copy()
            n_newton = it_r + it
        lam, dlam = 0.0, 1.0 / steps
        u = start.copy()
        while lam < 1.0 - 1e-12:
            lam_new = min(1.0, lam + dlam)
            trial = u.copy()
            # predictor: move the prescribed dofs, interior by the linearized response
            K = self.tangent(trial.reshape(-1, 3), P, alpha)
            _, R = self.energy_forces(trial.reshape(-1, 3), P, alpha)
            dup = (start + lam_new * (target - start) - trial)[fx]
            rhs = -(R[fr] + K[fr][:, fx] @ dup)
            trial[fx] += dup
            duf = spsolve(K[fr][:, fr].tocsc(), rhs)
            base = trial.copy()
            for sc in (1.0, 0.5, 0.25, 0.125, 0.0):                  # backtrack the predictor if it inverts elements
                trial[fr] = base[fr] + sc * duf
                if np.isfinite(self.energy_forces(trial.reshape(-1, 3), P, alpha)[0]):
                    break
            ok, it = self._newton(trial, P, alpha, tol, max_iter)
            n_newton += it
            if ok:
                u, lam = trial, lam_new
                n_steps += 1
                dlam = min(1.5 * dlam, 1.0 - lam) if it <= 5 else dlam
                if verbose:
                    print(f"    load {lam:.3f}: {it} Newton iterations ({time.time() - t0:.1f} s)", flush=True)
            else:
                dlam /= 2
                if dlam < 1e-3:
                    return u.reshape(-1, 3), dict(success=False, load=lam, newton=n_newton, steps=n_steps,
                                                  seconds=time.time() - t0)
        _, R = self.energy_forces(u.reshape(-1, 3), P, alpha)
        return u.reshape(-1, 3), dict(success=True, load=1.0, newton=n_newton, steps=n_steps,
                                      seconds=time.time() - t0, max_free_force=float(np.abs(R[fr]).max()),
                                      typical_reaction=float(np.median(np.abs(R[fx][np.abs(R[fx]) > 0])))
                                      if np.any(R[fx]) else 0.0)

    def _newton(self, u, P, alpha, tol, max_iter, dofs=None):
        """Newton iterations in place on u over the free dofs (or the subset dofs). Converged when the largest nodal
        force on those dofs is below tol times the largest reaction. Line search: backtracking on the energy (Armijo);
        when the energy change is at roundoff level, a step that lowers the residual is accepted."""
        fr, fx = (self.free if dofs is None else dofs), self.fixed
        E, R = self.energy_forces(u.reshape(-1, 3), P, alpha)
        if not np.isfinite(E):
            return False, 0
        target = tol * max(np.abs(R[fx]).max(), 1e-9)
        for it in range(1, max_iter + 1):
            if np.abs(R[fr]).max() < target:
                return True, it - 1
            K = self.tangent(u.reshape(-1, 3), P, alpha)
            du = spsolve(K[fr][:, fr].tocsc(), -R[fr])
            r0 = np.linalg.norm(R[fr])
            s = 1.0
            while s > 1e-4:
                trial = u.copy()
                trial[fr] += s * du
                E1, R1 = self.energy_forces(trial.reshape(-1, 3), P, alpha)
                if np.isfinite(E1):
                    if E1 <= E + 1e-4 * s * (R[fr] @ du):
                        break
                    if abs(E1 - E) <= 1e-11 * abs(E) and np.linalg.norm(R1[fr]) < r0:
                        break
                s /= 2
            if s <= 1e-4:
                return False, it
            u[:] = trial
            E, R = E1, R1
            target = tol * max(np.abs(R[fx]).max(), 1e-9)
        return bool(np.abs(R[fr]).max() < target), max_iter

    def from_membrane(self, xm):
        """Initial guess from a membrane solution xm (ny, nx, 2) on the same grid: the in-plane positions in every
        layer, the thickness scaled by the domain's mean in-plane area change about r = 0 (the top face or the mid-plane). Local thinning is
        left to solve(), which relaxes the through-thickness positions first (scaling each column by its own area
        change shears neighboring columns apart where the boundary is kinked)."""
        def cell_area(p):
            a = p[:-1, :-1]; b = p[:-1, 1:]; c = p[1:, 1:]; d = p[1:, :-1]
            return 0.5 * np.abs((c[..., 0] - a[..., 0]) * (d[..., 1] - b[..., 1])
                                - (d[..., 0] - b[..., 0]) * (c[..., 1] - a[..., 1]))
        xm = np.asarray(xm)
        J2 = cell_area(xm).sum() / cell_area(self.X[0, ..., :2]).sum()
        x = self.X.copy()
        x[..., :2] = xm[None]
        x[..., 2] = self.X[..., 2] / J2
        return x.reshape(-1, 3)

    def from_coarser(self, x_coarse, layers_coarse):
        """Initial guess from a solution with fewer element layers on the same grid and thickness: node positions
        interpolated linearly along each node column (mesh sequencing for the through-thickness refinement check)."""
        xc = np.asarray(x_coarse).reshape(layers_coarse + 1, self.ny, self.nx, 3)
        sc = np.linspace(0.0, 1.0, layers_coarse + 1)
        sf = np.linspace(0.0, 1.0, self.L + 1)
        k = np.clip(np.searchsorted(sc, sf, side="right") - 1, 0, layers_coarse - 1)
        w = ((sf - sc[k]) / (sc[k + 1] - sc[k]))[:, None, None, None]
        return ((1 - w) * xc[k] + w * xc[k + 1]).reshape(-1, 3)

    # ------------------------------------------------------------------ results
    def top(self, x):
        return x.reshape(self.X.shape)[-1]

    def cauchy(self, x, P, alpha=None):
        """Cauchy stress (kPa) at every Gauss point (ne, 8, 3, 3), including the element pressure."""
        alpha = np.ones(self.ne) if alpha is None else np.asarray(alpha, float).ravel()
        xe = x[self.conn]
        F = np.einsum("eai,egaj->egij", xe, self.gradN)
        J = np.linalg.det(F)
        theta = (J * self.detJ0).sum(1) / self.V0
        dU = (alpha * (-P[0]) + self._kappa(P) * np.log(theta)) / theta
        s = np.einsum("egij,egkj->egik", first_pk(F, P), F) / J[..., None, None] * alpha[:, None, None, None]
        return s + dU[:, None, None, None] * np.eye(3), F, theta


# ---------------------------------------------------------------------- self-checks
def _checks():
    import constitutive as M
    rng = np.random.default_rng(0)
    P = np.array([5.706, 4.072, 1.692, 9.883, 10.488, 0.991, 7.927, 45.573])     # N1739/test1 bulk fit
    # 1) first PK = dW/dF
    F = np.eye(3) + rng.normal(0, 0.1, (20, 3, 3)); F[:, 0, 0] += 0.25; F[:, 1, 1] += 0.12
    h = 1e-6
    num = np.empty_like(F)
    for i in range(3):
        for j in range(3):
            d = np.zeros((3, 3)); d[i, j] = h
            num[:, i, j] = (energy_density(F + d, P) - energy_density(F - d, P)) / (2 * h)
    print(f"1) first PK vs dW/dF: max rel error {np.abs(num - first_pk(F, P)).max() / np.abs(num).max():.1e}")

    # 2) element forces = d(energy)/dx on a small distorted slab
    Xg = np.stack(np.meshgrid(np.linspace(0, 2, 4), np.linspace(0, 1.6, 3)), -1)
    slab = Slab(Xg, 1.0 + 0.2 * rng.random((3, 4)), layers=2)
    x = slab.X.reshape(-1, 3) * [1.2, 1.1, 0.8] + rng.normal(0, 0.01, (slab.nn, 3))
    al = 1 + 0.3 * rng.random(slab.ne)
    E0, g = slab.energy_forces(x, P, al)
    num = np.empty(3 * slab.nn)
    for k in range(3 * slab.nn):
        xp, xm = x.ravel().copy(), x.ravel().copy()
        xp[k] += 1e-6; xm[k] -= 1e-6
        num[k] = (slab.energy_forces(xp.reshape(-1, 3), P, al)[0] - slab.energy_forces(xm.reshape(-1, 3), P, al)[0]) / 2e-6
    print(f"2) nodal forces vs d(energy)/dx: max rel error {np.abs(num - g).max() / np.abs(g).max():.1e}")

    # 3) homogeneous biaxial stretch: side faces stretched affinely, faces free -> uniform plane stress
    for lz, lt, shear in ((1.35, 1.17, 0.0), (1.30, 1.05, 0.0), (1.25, 1.12, 0.08)):
        Xg = np.stack(np.meshgrid(np.linspace(0, 4, 5), np.linspace(0, 3, 4)), -1)
        slab = Slab(Xg, 2.7, layers=2)
        Fin = np.array([[lz, shear], [0.0, lt]])
        xb = np.einsum("ij,...j->...i", Fin, Xg)
        x, info = slab.solve(xb, P, steps=4)
        s, F, theta = slab.cauchy(x, P)
        F2 = Fin
        t_m = M.cauchy(F2, P)                                         # membrane plane stress
        smean = s.reshape(-1, 3, 3).mean(0)
        spread = np.abs(s.reshape(-1, 3, 3) - smean).max()
        thick = (slab.top(x)[..., 2] - x.reshape(slab.X.shape)[0, ..., 2]).mean()
        print(f"3) F_inplane = [[{lz}, {shear}], [0, {lt}]]: solid t_zz {smean[0, 0]:.4f}, t_tt {smean[1, 1]:.4f}, "
              f"t_zt {smean[0, 1]:.4f}, t_rr {smean[2, 2]:.1e} kPa | membrane {t_m[0, 0]:.4f}, {t_m[1, 1]:.4f}, "
              f"{t_m[0, 1]:.4f} | spread {spread:.1e} | J-1 {np.abs(theta - 1).max():.1e} | thickness "
              f"{thick:.4f} vs {2.7 / (lz * lt):.4f} mm | {info['newton']} Newton it, {info['seconds']:.1f} s")


if __name__ == "__main__":
    _checks()
