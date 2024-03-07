#![feature(thread_id_value)]
use nsp_core::time_step;
use nsp_kerr_schild_kerr::KerrSchildKerr;

fn main() {
    let ksk = KerrSchildKerr::new(1.0, 0.5).unwrap();
    let h = time_step::compute_hamiltonian(ksk, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0).unwrap();
    println!("H = {}", h);
}
