#![feature(thread_id_value)]

use nalgebra::Matrix3;
use nalgebra::Vector3;
use nsp_core::error;
use nsp_core::log_error;
use nsp_core::metric::SpacetimeMetric;

#[derive(Debug)]
pub struct KerrSchildKerr {
    mass: f64,
    spin: f64,
}

impl KerrSchildKerr {
    pub fn new(mass: f64, spin: f64) -> Result<KerrSchildKerr, error::Error> {
        if mass * mass - spin * spin < 0.0 {
            log_error!(
                "Unphysical mass and spin combinations: M = {}, a = {}",
                mass,
                spin
            );
            Err(error::Error::UnphysicalMetric)
        } else {
            Ok(KerrSchildKerr { mass, spin })
        }
    }
}

fn r2(a: f64, x: f64, y: f64, z: f64) -> f64 {
    let rho2 = x * x + y * y + z * z;
    let rho2_m_a2 = rho2 - a * a;
    let a2_z2 = a * a * z * z;
    0.5 * rho2_m_a2 + (0.25 * rho2_m_a2 * rho2_m_a2 + a2_z2).sqrt()
}

fn h(a: f64, m: f64, x: f64, y: f64, z: f64) -> f64 {
    let r_2 = r2(a, x, y, z);
    let r = r_2.sqrt();
    let r3 = r * r_2;
    let r4 = r_2 * r_2;
    let a2z2 = a * a * z * z;
    m * r3 / (r4 + a2z2)
}

fn l_d(a: f64, x: f64, y: f64, z: f64) -> Vector3<f64> {
    let r_2 = r2(a, x, y, z);
    let r = r_2.sqrt();
    let r2_p_a2 = r_2 + a * a;
    Vector3::new((r * x + a * y) / r2_p_a2, (r * y - a * x) / r2_p_a2, z / r)
}

impl SpacetimeMetric for KerrSchildKerr {
    fn alpha(&self, _: f64, x: f64, y: f64, z: f64) -> f64 {
        let hh = h(self.spin, self.mass, x, y, z);
        let den = 1.0 / (1.0 + 2.0 * hh);
        den.sqrt()
    }

    fn beta_u(&self, _: f64, x: f64, y: f64, z: f64) -> Vector3<f64> {
        let hh = h(self.spin, self.mass, x, y, z);
        let den = 1.0 + 2.0 * hh;
        let prefactor = 2.0 * hh / den;
        prefactor * l_d(self.spin, x, y, z)
    }

    fn gamma_uu(&self, _: f64, x: f64, y: f64, z: f64) -> Matrix3<f64> {
        let hh = h(self.spin, self.mass, x, y, z);
        let den = 1.0 + 2.0 * hh;
        let prefactor = 2.0 * hh / den;
        let l = l_d(self.spin, x, y, z);

        let eta_ij = Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

        let h_ij = Matrix3::new(
            l[0] * l[0],
            l[0] * l[1],
            l[0] * l[2],
            l[1] * l[0],
            l[1] * l[1],
            l[1] * l[2],
            l[2] * l[0],
            l[2] * l[1],
            l[2] * l[2],
        );

        eta_ij - prefactor * h_ij
    }
}
