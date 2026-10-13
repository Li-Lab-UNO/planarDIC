"""
Bulk four-fiber family model for human femoropopliteal artery
(Jadidi et al. 2019, Biomech Model Mechanobiol 18:1591, Eqs. 1-3 and Appendix).

    W = Cgr/2 (I1 - 3)
      + sum_i  C1_i / (4 C2_i) * (exp(C2_i <I4_i - 1>^2) - 1),   i = el, smc, col1, col2

Fiber directions in the reference configuration (z = axial, theta = circumferential):
    el   : M = e_z
    smc  : M = e_theta
    col  : M = cos(g) e_z +/- sin(g) e_theta   (two symmetric families, same C1/C2)

Incompressible membrane under plane stress (t_rr = 0):
    t = -p I + Cgr B + sum_i C1_i <I4_i - 1> exp(C2_i (I4_i - 1)^2) m_i (x) m_i,   m_i = F M_i
    p = Cgr * lambda_r^2,  lambda_r = 1 / det(F_inplane)

Parameter vector (stresses in kPa, gamma in degrees):
    P = [Cgr, C1_el, C2_el, C1_smc, C2_smc, C1_col, C2_col, gamma_deg]
"""
import numpy as np

PARAM_NAMES = ("Cgr", "C1_el", "C2_el", "C1_smc", "C2_smc", "C1_col", "C2_col", "gamma_deg")
PARAM_UNITS = ("kPa", "kPa", "-", "kPa", "-", "kPa", "-", "deg")


def _fiber_dirs(gamma_deg):
    g = np.radians(gamma_deg)
    return [np.array([1.0, 0.0]), np.array([0.0, 1.0]),
            np.array([np.cos(g), np.sin(g)]), np.array([np.cos(g), -np.sin(g)])]


def _mac(x):
    return np.maximum(x, 0.0)


def strain_energy(F, P):
    """W (kPa) for in-plane deformation gradients F (..., 2, 2)."""
    Cgr, C1e, C2e, C1s, C2s, C1c, C2c, g = P
    C = np.einsum("...ki,...kj->...ij", F, F)
    lr2 = 1.0 / np.linalg.det(F) ** 2
    W = 0.5 * Cgr * (C[..., 0, 0] + C[..., 1, 1] + lr2 - 3)
    for M, c1, c2 in zip(_fiber_dirs(g), (C1e, C1s, C1c, C1c), (C2e, C2s, C2c, C2c)):
        I4 = np.einsum("i,...ij,j->...", M, C, M)
        W = W + c1 / (4 * c2) * (np.exp(c2 * _mac(I4 - 1) ** 2) - 1)
    return W


def cauchy(F, P):
    """
    In-plane Cauchy stress (kPa) of the incompressible membrane, plane stress.
    F: (..., 2, 2) in-plane deformation gradient, index 0 = z, 1 = theta.
    Returns t (..., 2, 2) with t[...,0,0] = t_zz, t[...,1,1] = t_thetatheta, t[...,0,1] = t_ztheta.
    """
    Cgr, C1e, C2e, C1s, C2s, C1c, C2c, g = P
    F = np.asarray(F, float)
    B = np.einsum("...ik,...jk->...ij", F, F)
    lr2 = 1.0 / np.linalg.det(F) ** 2
    t = Cgr * (B - lr2[..., None, None] * np.eye(2))       # isotropic part minus p = Cgr*lr^2
    for M, c1, c2 in zip(_fiber_dirs(g), (C1e, C1s, C1c, C1c), (C2e, C2s, C2c, C2c)):
        m = np.einsum("...ij,j->...i", F, M)
        e = np.einsum("...i,...i->...", m, m) - 1          # I4 - 1
        w = c1 * _mac(e) * np.exp(c2 * e ** 2)              # = 2 dW/dI4
        t = t + w[..., None, None] * np.einsum("...i,...j->...ij", m, m)
    return t


def cauchy_biaxial(lz, lt, P):
    """Closed form for F = diag(lz, lt) (Jadidi et al. 2019, Appendix). Returns t_zz, t_tt."""
    Cgr, C1e, C2e, C1s, C2s, C1c, C2c, g = P
    g = np.radians(g)
    lr2 = 1.0 / (lz * lt) ** 2
    I4c = lz ** 2 * np.cos(g) ** 2 + lt ** 2 * np.sin(g) ** 2
    col = 2 * C1c * _mac(I4c - 1) * np.exp(C2c * (I4c - 1) ** 2)
    tzz = Cgr * lz ** 2 + C1e * _mac(lz ** 2 - 1) * np.exp(C2e * (lz ** 2 - 1) ** 2) * lz ** 2 \
        + col * lz ** 2 * np.cos(g) ** 2 - Cgr * lr2
    ttt = Cgr * lt ** 2 + C1s * _mac(lt ** 2 - 1) * np.exp(C2s * (lt ** 2 - 1) ** 2) * lt ** 2 \
        + col * lt ** 2 * np.sin(g) ** 2 - Cgr * lr2
    return tzz, ttt


def first_pk_membrane(F, P, H):
    """
    Membrane tension per unit reference length, N = H t F^-T (kPa*mm = N/m), with H the
    reference thickness (mm). Current thickness h = H/J2 and 2D Nanson (J2 F^-T) give
    N = J2 (t h) F^-T = H t F^-T. Equilibrium on the reference grid: Div N = 0.
    """
    t = cauchy(F, P)
    return H * np.einsum("...ik,...jk->...ij", t, np.linalg.inv(F))


if __name__ == "__main__":
    # self-checks
    rng = np.random.default_rng(0)
    P = np.array([5.0, 40.0, 1.8, 10.9, 6.3, 5.4, 7.2, 50.6])   # Jadidi 2019 Table 1, 46.5 y group
    lz, lt = rng.uniform(1.0, 1.4, 200), rng.uniform(1.0, 1.3, 200)
    F = np.zeros((200, 2, 2)); F[:, 0, 0], F[:, 1, 1] = lz, lt
    t = cauchy(F, P); tzz, ttt = cauchy_biaxial(lz, lt, P)
    print("general vs closed form (diag F): max |diff| kPa",
          np.abs(t[:, 0, 0] - tzz).max(), np.abs(t[:, 1, 1] - ttt).max())
    # energy consistency, incompressible plane stress: t_zz = lz dW/dlz, t_tt = lt dW/dlt
    h = 1e-6
    def Wd(a, b):
        G = np.zeros(a.shape + (2, 2)); G[..., 0, 0], G[..., 1, 1] = a, b
        return strain_energy(G, P)
    dz = lz * (Wd(lz + h, lt) - Wd(lz - h, lt)) / (2 * h)
    dt = lt * (Wd(lz, lt + h) - Wd(lz, lt - h)) / (2 * h)
    print("stress vs dW/dlambda: max rel diff", np.max(np.abs(dz - tzz) / np.abs(tzz)),
          np.max(np.abs(dt - ttt) / np.abs(ttt)))
    # general F with shear: power-conjugacy check  dW = t : (dF F^-1) for incompressible membrane
    A = np.eye(2) + rng.normal(0, .1, (200, 2, 2)); A[:, 0, 0] += .2
    dF = rng.normal(0, 1, (200, 2, 2)) * 1e-6
    dW = strain_energy(A + dF, P) - strain_energy(A - dF, P)
    L = np.einsum("nij,njk->nik", 2 * dF, np.linalg.inv(A))
    pred = np.einsum("nij,nij->n", cauchy(A, P), L)      # t : sym(L) (t symmetric)
    print("general F (with shear) energy consistency: max rel diff",
          np.max(np.abs(dW - pred) / (np.abs(pred) + 1e-12)))
