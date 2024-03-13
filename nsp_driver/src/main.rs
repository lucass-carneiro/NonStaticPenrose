#![feature(thread_id_value)]
use nsp_core::{derivatives, time_step};
use nsp_kerr_schild_kerr::KerrSchildKerr;

fn main() {
    let ksk = KerrSchildKerr::new(1.0, 0.5).unwrap();

    let s = nsp_core::time_step::State {
        px: 1.0,
        py: 1.0,
        pz: 1.0,
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };

    let h = time_step::compute_hamiltonian(&ksk, &s, 1.0, 0.0).unwrap();
    let hp = derivatives::fd_4(
        &ksk,
        &s,
        1.0,
        0.0,
        derivatives::DerivativeComponent::X,
        1.0e-3,
    )
    .unwrap();

    println!("H = {}", h);
    println!("dH/dpx = {}", hp);
}
