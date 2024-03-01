use nalgebra::Matrix3;
use nalgebra::Vector3;
use nsp_core::metric::SpacetimeMetric;

#[derive(Debug)]
pub struct KerrSchildSchwarzschild {
    mass: f64,
}

impl KerrSchildSchwarzschild {
    pub fn new(mass: f64) -> KerrSchildSchwarzschild {
        KerrSchildSchwarzschild { mass }
    }
}

fn l_d(x: f64, y: f64, z: f64) -> Vector3<f64> {
    let r = (x * x + y * y + z * z).sqrt();
    Vector3::new(x / r, y / r, z / r)
}

impl SpacetimeMetric for KerrSchildSchwarzschild {
    fn alpha(&self, _: f64, x: f64, y: f64, z: f64) -> f64 {
        let r = (x * x + y * y + z * z).sqrt();
        (1.0 / (1.0 + 2.0 * self.mass / r)).sqrt()
    }

    fn beta_u(&self, t: f64, x: f64, y: f64, z: f64) -> Vector3<f64> {
        let r = (x * x + y * y + z * z).sqrt();
        let alp = self.alpha(t, x, y, z);
        let l = l_d(x, y, z);
        let prefactor = 2.0 * self.mass / r * alp * alp;
        prefactor * l
    }

    fn gamma_dd(&self, _: f64, x: f64, y: f64, z: f64) -> Matrix3<f64> {
        let r = (x * x + y * y + z * z).sqrt();
        let l = l_d(x, y, z);

        let eta_ij = Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

        let pf = 2.0 * self.mass / r;

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

        eta_ij - pf * h_ij
    }

    fn k_dd(&self, t: f64, x: f64, y: f64, z: f64) -> Matrix3<f64> {
        let r2 = x * x + y * y + z * z;
        let r = r2.sqrt();

        let l = l_d(x, y, z);
        let a = self.alpha(t, x, y, z);

        let pf = 2.0 * self.mass * a / r2;
        let pf2 = 2.0 + self.mass / r;

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

        pf * (eta_ij - pf2 * h_ij)
    }

    fn k(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        let r2 = x * x + y * y + z * z;
        let r = r2.sqrt();

        let a = self.alpha(t, x, y, z);
        let a3 = a * a * a;

        (2.0 * self.mass * a3 / r2) * (1.0 + 3.0 * self.mass / r)
    }
}
