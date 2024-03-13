use crate::{
    metric::SpacetimeMetric,
    time_step::{compute_hamiltonian, State},
};

pub enum DerivativeComponent {
    Px,
    Py,
    Pz,
    X,
    Y,
    Z,
}

pub fn fd_2<T: SpacetimeMetric>(
    sm: &T,
    state: &State,
    m: f64,
    t: f64,
    c: DerivativeComponent,
    h: f64,
) -> Result<f64, crate::error::Error> {
    let mut state_p = state.clone();
    let mut state_m = state.clone();

    match c {
        DerivativeComponent::Px => {
            state_p.px += h;
            state_m.px -= h;
        }

        DerivativeComponent::Py => {
            state_p.py += h;
            state_p.py -= h;
        }

        DerivativeComponent::Pz => {
            state_p.pz += h;
            state_m.pz -= h;
        }

        DerivativeComponent::X => {
            state_p.x += h;
            state_m.x -= h;
        }

        DerivativeComponent::Y => {
            state_p.y += h;
            state_m.y -= h;
        }

        DerivativeComponent::Z => {
            state_p.z += h;
            state_m.z -= h;
        }
    }

    let h_p = compute_hamiltonian(sm, &state_p, m, t)?;
    let h_m = compute_hamiltonian(sm, &state_m, m, t)?;

    Ok((h_p - h_m) / (2.0 * h))
}

pub fn fd_4<T: SpacetimeMetric>(
    sm: &T,
    state: &State,
    m: f64,
    t: f64,
    c: DerivativeComponent,
    h: f64,
) -> Result<f64, crate::error::Error> {
    let mut state_p = state.clone();
    let mut state_pp = state.clone();
    let mut state_m = state.clone();
    let mut state_mm = state.clone();

    match c {
        DerivativeComponent::Px => {
            state_p.px += h;
            state_pp.px += 2.0 * h;
            state_m.px -= h;
            state_mm.px -= 2.0 * h;
        }

        DerivativeComponent::Py => {
            state_p.py += h;
            state_pp.py += 2.0 * h;
            state_m.py -= h;
            state_mm.py -= 2.0 * h;
        }

        DerivativeComponent::Pz => {
            state_p.pz += h;
            state_pp.pz += 2.0 * h;
            state_m.pz -= h;
            state_mm.pz -= 2.0 * h;
        }

        DerivativeComponent::X => {
            state_p.x += h;
            state_pp.x += 2.0 * h;
            state_m.x -= h;
            state_mm.x -= 2.0 * h;
        }

        DerivativeComponent::Y => {
            state_p.y += h;
            state_pp.y += 2.0 * h;
            state_m.y -= h;
            state_mm.y -= 2.0 * h;
        }

        DerivativeComponent::Z => {
            state_p.z += h;
            state_pp.z += 2.0 * h;
            state_m.z -= h;
            state_mm.z -= 2.0 * h;
        }
    }

    let h_p = compute_hamiltonian(sm, &state_p, m, t)?;
    let h_pp = compute_hamiltonian(sm, &state_pp, m, t)?;
    let h_m = compute_hamiltonian(sm, &state_m, m, t)?;
    let h_mm = compute_hamiltonian(sm, &state_mm, m, t)?;

    Ok((h_mm - 8.0 * h_m + 8.0 * h_p - h_pp) / (12.0 * h))
}
