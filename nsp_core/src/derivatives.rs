use crate::{
    metric::SpacetimeMetric,
    time_step::{compute_hamiltonian, State},
};

enum DerivativeComponent {
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
    //c: DerivativeComponent,
    h: f64,
) -> Result<f64, crate::error::Error> {
    let state_p = State {
        px: state.px + h,
        py: state.py,
        pz: state.pz,
        x: state.x,
        y: state.y,
        z: state.z,
    };

    let state_m = State {
        px: state.px - h,
        py: state.py,
        pz: state.pz,
        x: state.x,
        y: state.y,
        z: state.z,
    };

    let h = compute_hamiltonian(sm, state, m, t)?;
    let h_p = compute_hamiltonian(sm, &state_p, m, t)?;
    let h_m = compute_hamiltonian(sm, &state_m, m, t)?;

    Ok((h_p - h_m) / (2.0 * h))
}
