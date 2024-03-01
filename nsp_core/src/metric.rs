use nalgebra::Matrix3;
use nalgebra::Vector3;

const FD_STEP: f64 = 1.0e-2;

const FD_8_WEIGTHS: [f64; 9] = [
    1.0 / (280.0 * FD_STEP),
    -4.0 / (105.0 * FD_STEP),
    1.0 / (5.0 * FD_STEP),
    -4.0 / (5.0 * FD_STEP),
    0.0,
    4.0 / (5.0 * FD_STEP),
    -1.0 / (5.0 * FD_STEP),
    4.0 / (105.0 * FD_STEP),
    -1.0 / (280.0 * FD_STEP),
];

pub trait SpacetimeMetric {
    fn alpha(&self, t: f64, x: f64, y: f64, z: f64) -> f64;

    fn beta_u(&self, t: f64, x: f64, y: f64, z: f64) -> Vector3<f64>;

    fn gamma_dd(&self, t: f64, x: f64, y: f64, z: f64) -> Matrix3<f64>;

    fn k_dd(&self, t: f64, x: f64, y: f64, z: f64) -> Matrix3<f64>;
    fn k(&self, t: f64, x: f64, y: f64, z: f64) -> f64;

    fn d_alpha_dt(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        let fn_values = [
            self.alpha(t - 4.0 * FD_STEP, x, y, z),
            self.alpha(t - 3.0 * FD_STEP, x, y, z),
            self.alpha(t - 2.0 * FD_STEP, x, y, z),
            self.alpha(t - 1.0 * FD_STEP, x, y, z),
            self.alpha(t - 0.0 * FD_STEP, x, y, z),
            self.alpha(t + 1.0 * FD_STEP, x, y, z),
            self.alpha(t + 2.0 * FD_STEP, x, y, z),
            self.alpha(t + 3.0 * FD_STEP, x, y, z),
            self.alpha(t + 4.0 * FD_STEP, x, y, z),
        ];

        let mut fd = 0.0;

        for i in 0..9 {
            fd += FD_8_WEIGTHS[i] * fn_values[i];
        }

        fd
    }

    fn d_alpha_dx(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        let fn_values = [
            self.alpha(t, x - 4.0 * FD_STEP, y, z),
            self.alpha(t, x - 3.0 * FD_STEP, y, z),
            self.alpha(t, x - 2.0 * FD_STEP, y, z),
            self.alpha(t, x - 1.0 * FD_STEP, y, z),
            self.alpha(t, x - 0.0 * FD_STEP, y, z),
            self.alpha(t, x + 1.0 * FD_STEP, y, z),
            self.alpha(t, x + 2.0 * FD_STEP, y, z),
            self.alpha(t, x + 3.0 * FD_STEP, y, z),
            self.alpha(t, x + 4.0 * FD_STEP, y, z),
        ];

        let mut fd = 0.0;

        for i in 0..9 {
            fd += FD_8_WEIGTHS[i] * fn_values[i];
        }

        fd
    }

    fn d_alpha_dy(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        let fn_values = [
            self.alpha(t, x, y - 4.0 * FD_STEP, z),
            self.alpha(t, x, y - 3.0 * FD_STEP, z),
            self.alpha(t, x, y - 2.0 * FD_STEP, z),
            self.alpha(t, x, y - 1.0 * FD_STEP, z),
            self.alpha(t, x, y - 0.0 * FD_STEP, z),
            self.alpha(t, x, y + 1.0 * FD_STEP, z),
            self.alpha(t, x, y + 2.0 * FD_STEP, z),
            self.alpha(t, x, y + 3.0 * FD_STEP, z),
            self.alpha(t, x, y + 4.0 * FD_STEP, z),
        ];

        let mut fd = 0.0;

        for i in 0..9 {
            fd += FD_8_WEIGTHS[i] * fn_values[i];
        }

        fd
    }

    fn d_alpha_dz(&self, t: f64, x: f64, y: f64, z: f64) -> f64 {
        let fn_values = [
            self.alpha(t, x, y, z - 4.0 * FD_STEP),
            self.alpha(t, x, y, z - 3.0 * FD_STEP),
            self.alpha(t, x, y, z - 2.0 * FD_STEP),
            self.alpha(t, x, y, z - 1.0 * FD_STEP),
            self.alpha(t, x, y, z - 0.0 * FD_STEP),
            self.alpha(t, x, y, z + 1.0 * FD_STEP),
            self.alpha(t, x, y, z + 2.0 * FD_STEP),
            self.alpha(t, x, y, z + 3.0 * FD_STEP),
            self.alpha(t, x, y, z + 4.0 * FD_STEP),
        ];

        let mut fd = 0.0;

        for i in 0..9 {
            fd += FD_8_WEIGTHS[i] * fn_values[i];
        }

        fd
    }
}
