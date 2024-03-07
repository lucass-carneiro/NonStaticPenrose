use nalgebra::Matrix3;
use nalgebra::Vector3;

pub trait SpacetimeMetric {
    fn alpha(&self, t: f64, x: f64, y: f64, z: f64) -> f64;
    fn beta_u(&self, t: f64, x: f64, y: f64, z: f64) -> Vector3<f64>;
    fn gamma_uu(&self, t: f64, x: f64, y: f64, z: f64) -> Matrix3<f64>;
}
