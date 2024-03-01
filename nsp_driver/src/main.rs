#![feature(thread_id_value)]
use nsp_core::metric::SpacetimeMetric;
use nsp_kerr_schild_schwarzschild::KerrSchildSchwarzschild;

fn main() {
    let is = KerrSchildSchwarzschild::new(1.0);
    nsp_core::log_info!("Hey! {}", is.alpha(0.0, 1.0, 1.0, 1.0));
    nsp_core::log_info!("Hey! {}", is.d_alpha_dx(0.0, 1.0, 1.0, 1.0));
    nsp_core::log_info!("Hey! {}", is.d_alpha_dy(0.0, 1.0, 1.0, 1.0));
    nsp_core::log_info!("Hey! {}", is.d_alpha_dz(0.0, 1.0, 1.0, 1.0));
}
