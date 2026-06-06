//! Repository orbits around vessel center.
//!
//! Repos orbit their vessel following Kepler-like dynamics:
//! core repos (small orbital radius) orbit faster, peripheral
//! repos orbit slower. Orbital period obeys Kepler's third law:
//! T² ∝ r³.

use serde::{Deserialize, Serialize};

/// A repository orbiting its vessel center.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Repo {
    pub name: String,
    pub vessel: String,
    /// Orbital radius from vessel center in dependency-space.
    pub orbital_radius: f64,
    /// Current angle (radians).
    pub angle: f64,
    /// Angular velocity (radians per time step).
    pub angular_velocity: f64,
}

impl Repo {
    /// Create a new repo at the given orbital radius.
    /// Angular velocity is computed from Kepler's 3rd law analog:
    /// ω = √(μ / r³) where μ = G * vessel_mass.
    pub fn new(name: impl Into<String>, vessel: impl Into<String>, orbital_radius: f64, mu: f64) -> Self {
        let angular_velocity = if orbital_radius > 0.0 && mu > 0.0 {
            (mu / orbital_radius.powi(3)).sqrt()
        } else {
            0.0
        };
        Self {
            name: name.into(),
            vessel: vessel.into(),
            orbital_radius,
            angle: 0.0,
            angular_velocity,
        }
    }

    /// Create a repo with explicit angular velocity (override Kepler).
    pub fn with_angular_velocity(name: impl Into<String>, vessel: impl Into<String>, orbital_radius: f64, angle: f64, angular_velocity: f64) -> Self {
        Self { name: name.into(), vessel: vessel.into(), orbital_radius, angle, angular_velocity }
    }

    /// Orbital period: T = 2π / ω
    pub fn period(&self) -> f64 {
        if self.angular_velocity.abs() < 1e-15 {
            f64::INFINITY
        } else {
            2.0 * std::f64::consts::PI / self.angular_velocity
        }
    }

    /// Advance the repo by one time step dt.
    pub fn step(&mut self, dt: f64) {
        self.angle += self.angular_velocity * dt;
        // Keep angle in [0, 2π)
        self.angle = self.angle.rem_euclid(2.0 * std::f64::consts::PI);
    }

    /// Cartesian position relative to vessel center.
    pub fn position(&self) -> (f64, f64) {
        let x = self.orbital_radius * self.angle.cos();
        let y = self.orbital_radius * self.angle.sin();
        (x, y)
    }

    /// Orbital speed: v = ω * r
    pub fn orbital_speed(&self) -> f64 {
        self.angular_velocity * self.orbital_radius
    }

    /// Orbital energy: E = -μ/(2r) (bound orbit)
    pub fn orbital_energy(&self, mu: f64) -> f64 {
        if self.orbital_radius > 0.0 {
            -mu / (2.0 * self.orbital_radius)
        } else {
            0.0
        }
    }

    /// Verify Kepler's 3rd law: T² * μ = 4π² * r³
    pub fn verify_kepler(&self, mu: f64) -> bool {
        if self.orbital_radius <= 0.0 || mu <= 0.0 { return true; }
        let t = self.period();
        let lhs = t * t * mu;
        let rhs = 4.0 * std::f64::consts::PI.powi(2) * self.orbital_radius.powi(3);
        (lhs - rhs).abs() < 1e-6 * rhs.abs().max(1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repo_orbital_period_kepler() {
        // T² = 4π²r³/μ → ω = √(μ/r³) → T = 2π/ω = 2π√(r³/μ)
        let r = Repo::new("core-lib", "Forgemaster", 2.0, 50.0);
        let expected_t = 2.0 * std::f64::consts::PI * (2.0_f64.powi(3) / 50.0).sqrt();
        assert!((r.period() - expected_t).abs() < 1e-10);
    }

    #[test]
    fn kepler_law_verification() {
        let r = Repo::new("test-repo", "CCC", 5.0, 100.0);
        assert!(r.verify_kepler(100.0));
    }

    #[test]
    fn repo_step_advances_angle() {
        let mut r = Repo::new("test", "V", 1.0, 1.0);
        let initial = r.angle;
        r.step(0.1);
        assert!(r.angle > initial);
    }

    #[test]
    fn repo_position_cartesian() {
        let r = Repo::with_angular_velocity("test", "V", 2.0, 0.0, 1.0);
        let (x, y) = r.position();
        assert!((x - 2.0).abs() < 1e-10);
        assert!(y.abs() < 1e-10);
    }

    #[test]
    fn repo_position_at_pi_over_2() {
        let r = Repo::with_angular_velocity("test", "V", 1.0, std::f64::consts::FRAC_PI_2, 1.0);
        let (x, y) = r.position();
        assert!(x.abs() < 1e-10);
        assert!((y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn orbital_speed() {
        let r = Repo::with_angular_velocity("test", "V", 3.0, 0.0, 2.0);
        assert!((r.orbital_speed() - 6.0).abs() < 1e-10);
    }

    #[test]
    fn closer_repos_orbit_faster() {
        let mu = 100.0;
        let inner = Repo::new("inner", "V", 1.0, mu);
        let outer = Repo::new("outer", "V", 10.0, mu);
        assert!(inner.angular_velocity > outer.angular_velocity);
    }

    #[test]
    fn serde_roundtrip_repo() {
        let r = Repo::with_angular_velocity("core", "Forgemaster", 3.5, 1.2, 0.8);
        let json = serde_json::to_string(&r).unwrap();
        let r2: Repo = serde_json::from_str(&json).unwrap();
        assert_eq!(r, r2);
    }
}
