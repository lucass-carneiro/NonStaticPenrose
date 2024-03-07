use crate::log_error;

use crate::metric::SpacetimeMetric;

pub fn compute_hamiltonian<T: SpacetimeMetric>(
    sm: T,
    m: f64,
    px: f64,
    py: f64,
    pz: f64,
    t: f64,
    x: f64,
    y: f64,
    z: f64,
) -> Result<f64, crate::error::Error> {
    let alpha = sm.alpha(t, x, y, z);
    let a2 = alpha * alpha;
    let beta_u = sm.beta_u(t, x, y, z);
    let g_uu = sm.gamma_uu(t, x, y, z);

    // gamma^{i j} p_i p_j
    let p1 = g_uu[(0, 0)] * px * px
        + 2.0 * g_uu[(0, 1)] * px * py
        + 2.0 * g_uu[(0, 2)] * px * pz
        + g_uu[(1, 1)] * py * py
        + 2.0 * g_uu[(1, 2)] * py * pz
        + g_uu[(2, 2)] * pz * pz;

    // beta^i beta^j p_i p_j
    let p2 = beta_u[0] * beta_u[0] * px * px
        + 2.0 * beta_u[0] * beta_u[1] * px * py
        + 2.0 * beta_u[0] * beta_u[2] * px * pz
        + beta_u[1] * beta_u[1] * py * py
        + 2.0 * beta_u[1] * beta_u[2] * py * pz
        + beta_u[2] * beta_u[2] * pz * pz;

    // beta^i p_i
    let p3 = beta_u[0] * px + beta_u[1] * py + beta_u[2] * pz;

    let p4 = a2 * (p1 + m * m) + p3 * p3 - p2;

    if p4 < 0.0 {
        log_error!(
            "Unphysical Hamiltonian: x = ({}, {}, {}, {}) p = ({}, {}, {})",
            t,
            x,
            y,
            z,
            px,
            py,
            pz
        );
        Err(crate::error::Error::UnphysicalHamiltonian)
    } else {
        Ok(p4.sqrt() - p3)
    }
}
