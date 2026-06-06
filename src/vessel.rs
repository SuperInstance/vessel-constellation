//! A vessel as a gravitational body.
//!
//! Each vessel has mass proportional to its repository count,
//! a position in dependency-space, and a velocity representing
//! its growth rate.

use serde::{Deserialize, Serialize};

/// A vessel is a gravitational body in the constellation.
///
/// Mass = number of repositories. Position lives in a
/// dependency-space where dimensions represent categories
/// of inter-project coupling.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Vessel {
    pub name: String,
    pub mass: f64,
    /// Position in dependency-space (dimensions = coupling categories).
    pub position: Vec<f64>,
    /// Velocity (growth rate vector).
    pub velocity: Vec<f64>,
}

impl Vessel {
    /// Create a new vessel at the origin with zero velocity.
    pub fn new(name: impl Into<String>, mass: f64, dimensions: usize) -> Self {
        Self {
            name: name.into(),
            mass,
            position: vec![0.0; dimensions],
            velocity: vec![0.0; dimensions],
        }
    }

    /// Create a vessel at a specific position.
    pub fn at(name: impl Into<String>, mass: f64, position: Vec<f64>) -> Self {
        let velocity = vec![0.0; position.len()];
        Self { name: name.into(), mass, position, velocity }
    }

    /// Create a vessel with position and velocity.
    pub fn with_velocity(name: impl Into<String>, mass: f64, position: Vec<f64>, velocity: Vec<f64>) -> Self {
        assert_eq!(position.len(), velocity.len(), "position and velocity must have same dimensionality");
        Self { name: name.into(), mass, position, velocity }
    }

    /// Number of spatial dimensions this vessel lives in.
    pub fn dimensions(&self) -> usize {
        self.position.len()
    }

    /// Euclidean distance to another vessel.
    pub fn distance_to(&self, other: &Vessel) -> f64 {
        self.position.iter()
            .zip(other.position.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Kinetic energy: KE = ½mv²
    pub fn kinetic_energy(&self) -> f64 {
        let v_sq: f64 = self.velocity.iter().map(|v| v * v).sum();
        0.5 * self.mass * v_sq
    }

    /// Momentum: p = mv
    pub fn momentum(&self) -> Vec<f64> {
        self.velocity.iter().map(|v| self.mass * v).collect()
    }

    /// Angular momentum about the origin: L = m(r × v).
    /// For 2D, returns a 1-element vec with the z-component.
    /// For 3D, returns the full cross product.
    pub fn angular_momentum(&self) -> Vec<f64> {
        cross_product(&self.position, &self.velocity)
            .iter()
            .map(|c| self.mass * c)
            .collect()
    }

    /// Center of mass between two vessels.
    pub fn center_of_mass(&self, other: &Vessel) -> Vec<f64> {
        let total = self.mass + other.mass;
        self.position.iter()
            .zip(other.position.iter())
            .map(|(a, b)| (self.mass * a + other.mass * b) / total)
            .collect()
    }
}

/// Cross product generalized: always returns 3-component vector.
/// If inputs are 2D, treats them as (x, y, 0).
pub fn cross_product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (ax, ay, az) = match a.len() {
        2 => (a[0], a[1], 0.0),
        3 => (a[0], a[1], a[2]),
        _ => (a.first().copied().unwrap_or(0.0), a.get(1).copied().unwrap_or(0.0), a.get(2).copied().unwrap_or(0.0)),
    };
    let (bx, by, bz) = match b.len() {
        2 => (b[0], b[1], 0.0),
        3 => (b[0], b[1], b[2]),
        _ => (b.first().copied().unwrap_or(0.0), b.get(1).copied().unwrap_or(0.0), b.get(2).copied().unwrap_or(0.0)),
    };
    vec![ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vessel_new_has_zero_position_and_velocity() {
        let v = Vessel::new("Forgemaster", 330.0, 3);
        assert_eq!(v.position, vec![0.0, 0.0, 0.0]);
        assert_eq!(v.velocity, vec![0.0, 0.0, 0.0]);
        assert_eq!(v.mass, 330.0);
    }

    #[test]
    fn vessel_mass_from_repo_count() {
        let fm = Vessel::new("Forgemaster", 330.0, 2);
        let ccc = Vessel::new("CCC", 116.0, 2);
        let jc = Vessel::new("JetsonClaw1", 76.0, 2);
        let oracle = Vessel::new("Oracle", 43.0, 2);
        assert_eq!(fm.mass, 330.0);
        assert_eq!(ccc.mass, 116.0);
        assert_eq!(jc.mass, 76.0);
        assert_eq!(oracle.mass, 43.0);
    }

    #[test]
    fn distance_between_vessels() {
        let a = Vessel::at("A", 1.0, vec![0.0, 0.0]);
        let b = Vessel::at("B", 1.0, vec![3.0, 4.0]);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn kinetic_energy_stationary() {
        let v = Vessel::at("V", 10.0, vec![1.0, 2.0]);
        assert!((v.kinetic_energy() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn kinetic_energy_moving() {
        let v = Vessel::with_velocity("V", 2.0, vec![0.0, 0.0], vec![3.0, 4.0]);
        // KE = 0.5 * 2.0 * (9 + 16) = 25.0
        assert!((v.kinetic_energy() - 25.0).abs() < 1e-10);
    }

    #[test]
    fn momentum_calculation() {
        let v = Vessel::with_velocity("V", 3.0, vec![0.0, 0.0], vec![2.0, 5.0]);
        assert_eq!(v.momentum(), vec![6.0, 15.0]);
    }

    #[test]
    fn center_of_mass_equal_mass() {
        let a = Vessel::at("A", 1.0, vec![0.0, 0.0]);
        let b = Vessel::at("B", 1.0, vec![2.0, 0.0]);
        let com = a.center_of_mass(&b);
        assert!((com[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn serde_roundtrip_vessel() {
        let v = Vessel::with_velocity("Forgemaster", 330.0, vec![1.0, 2.0, 3.0], vec![0.1, 0.2, 0.3]);
        let json = serde_json::to_string(&v).unwrap();
        let v2: Vessel = serde_json::from_str(&json).unwrap();
        assert_eq!(v, v2);
    }
}
