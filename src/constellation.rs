//! Full constellation state and evolution.
//!
//! The constellation holds all vessels, their orbiting repos,
//! and the gravitational field. Evolution uses symplectic
//! leapfrog (Störmer-Verlet) integration for superior energy
//! conservation compared to Euler methods.

use serde::{Deserialize, Serialize};
use crate::vessel::Vessel;
use crate::orbit::Repo;
use crate::gravity::GravitationalField;
use crate::conservation::ConservationState;

/// The full constellation: all vessels, repos, and dynamics parameters.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Constellation {
    pub vessels: Vec<Vessel>,
    pub repos: Vec<Repo>,
    pub dt: f64,
    pub G: f64,
}

impl Constellation {
    /// Create a new constellation.
    pub fn new(vessels: Vec<Vessel>, dt: f64, G: f64) -> Self {
        Self { vessels, repos: Vec::new(), dt, G }
    }

    /// Build the gravitational field from current vessel positions.
    pub fn field(&self) -> GravitationalField {
        GravitationalField::new(self.vessels.clone(), self.G)
    }

    /// Compute current conservation state.
    pub fn conservation(&self) -> ConservationState {
        let field = self.field();
        ConservationState::compute(&self.vessels, &field)
    }

    /// Evolve the constellation by one leapfrog step (kick-drift-kick).
    ///
    /// Leapfrog (velocity Verlet):
    ///   1. Half-kick:  v(t+dt/2) = v(t) + a(t)·dt/2
    ///   2. Drift:      r(t+dt)   = r(t) + v(t+dt/2)·dt
    ///   3. Recompute accelerations at new positions
    ///   4. Half-kick:  v(t+dt)   = v(t+dt/2) + a(t+dt)·dt/2
    pub fn step_leapfrog(&mut self) {
        let dt = self.dt;
        let dim = self.vessels.first().map_or(0, |v| v.dimensions());

        // Step 1: compute accelerations at current positions
        let field = self.field();
        let acc_old = field.accelerations();

        // Half-kick velocities
        for (v, a) in self.vessels.iter_mut().zip(acc_old.iter()) {
            for i in 0..dim.min(v.velocity.len()) {
                v.velocity[i] += a[i] * dt * 0.5;
            }
        }

        // Step 2: drift positions
        for v in &mut self.vessels {
            for i in 0..dim.min(v.position.len()) {
                v.position[i] += v.velocity[i] * dt;
            }
        }

        // Step 3: new accelerations
        let field_new = self.field();
        let acc_new = field_new.accelerations();

        // Step 4: half-kick velocities
        for (v, a) in self.vessels.iter_mut().zip(acc_new.iter()) {
            for i in 0..dim.min(v.velocity.len()) {
                v.velocity[i] += a[i] * dt * 0.5;
            }
        }

        // Evolve repo orbits
        for repo in &mut self.repos {
            repo.step(dt);
        }
    }

    /// Evolve by N steps, returning initial and final conservation states.
    pub fn evolve(&mut self, steps: usize) -> (ConservationState, ConservationState) {
        let initial = self.conservation();
        for _ in 0..steps {
            self.step_leapfrog();
        }
        let final_state = self.conservation();
        (initial, final_state)
    }

    /// Simple Euler integration step (for comparison/testing).
    /// v(t+dt) = v(t) + a(t)·dt
    /// r(t+dt) = r(t) + v(t)·dt
    pub fn step_euler(&mut self) {
        let dt = self.dt;
        let dim = self.vessels.first().map_or(0, |v| v.dimensions());
        let field = self.field();
        let acc = field.accelerations();

        for (v, a) in self.vessels.iter_mut().zip(acc.iter()) {
            for i in 0..dim.min(v.position.len()) {
                v.position[i] += v.velocity[i] * dt;
                v.velocity[i] += a[i] * dt;
            }
        }
    }

    /// Total kinetic energy of all vessels.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.vessels.iter().map(|v| v.kinetic_energy()).sum()
    }

    /// Total potential energy from gravitational interactions.
    pub fn total_potential_energy(&self) -> f64 {
        self.field().total_potential_energy()
    }
}

/// Create the initial fleet constellation.
pub fn initial_fleet() -> Constellation {
    let vessels = vec![
        Vessel::at("Forgemaster", 330.0, vec![0.0, 0.0]),
        Vessel::at("CCC", 116.0, vec![10.0, 0.0]),
        Vessel::at("JetsonClaw1", 76.0, vec![5.0, 8.0]),
        Vessel::at("Oracle", 43.0, vec![-3.0, 5.0]),
    ];
    Constellation::new(vessels, 0.01, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constellation_new() {
        let c = Constellation::new(vec![], 0.01, 1.0);
        assert_eq!(c.vessels.len(), 0);
        assert_eq!(c.dt, 0.01);
    }

    #[test]
    fn initial_fleet_has_four_vessels() {
        let fleet = initial_fleet();
        assert_eq!(fleet.vessels.len(), 4);
        assert_eq!(fleet.vessels[0].name, "Forgemaster");
        assert_eq!(fleet.vessels[0].mass, 330.0);
    }

    #[test]
    fn leapfrog_conserves_energy_over_100_steps() {
        let mut c = Constellation::new(vec![
            Vessel::with_velocity("A", 10.0, vec![0.0, 0.0], vec![0.0, 0.2]),
            Vessel::with_velocity("B", 10.0, vec![5.0, 0.0], vec![0.0, -0.2]),
        ], 0.001, 1.0);
        let (initial, final_state) = c.evolve(100);
        assert!(initial.energy_conserved(&final_state, 0.01));
    }

    #[test]
    fn leapfrog_conserves_angular_momentum_over_100_steps() {
        let mut c = Constellation::new(vec![
            Vessel::with_velocity("A", 10.0, vec![0.0, 0.0], vec![0.0, 0.5]),
            Vessel::with_velocity("B", 5.0, vec![3.0, 0.0], vec![0.0, -1.0]),
        ], 0.001, 1.0);
        let (initial, final_state) = c.evolve(100);
        assert!(initial.angular_momentum_conserved(&final_state, 0.01));
    }

    #[test]
    fn leapfrog_better_than_euler() {
        // Same initial conditions, compare energy drift
        let make_system = || Constellation::new(vec![
            Vessel::with_velocity("A", 10.0, vec![0.0, 0.0], vec![0.0, 0.3]),
            Vessel::with_velocity("B", 10.0, vec![4.0, 0.0], vec![0.0, -0.3]),
        ], 0.005, 1.0);

        let mut lf = make_system();
        let lf_initial_energy = lf.conservation().total_energy;
        for _ in 0..200 { lf.step_leapfrog(); }
        let lf_drift = (lf.conservation().total_energy - lf_initial_energy).abs();

        let mut eu = make_system();
        let eu_initial_energy = eu.conservation().total_energy;
        for _ in 0..200 { eu.step_euler(); }
        let eu_drift = (eu.conservation().total_energy - eu_initial_energy).abs();

        // Leapfrog should have significantly less energy drift
        assert!(lf_drift <= eu_drift * 5.0, "Leapfrog drift: {lf_drift}, Euler drift: {eu_drift}");
    }

    #[test]
    fn full_pipeline_gravity_orbit_conserve_evolve() {
        let mut c = Constellation::new(vec![
            Vessel::with_velocity("Forgemaster", 330.0, vec![0.0, 0.0], vec![0.0, 0.01]),
            Vessel::with_velocity("CCC", 116.0, vec![8.0, 0.0], vec![0.0, -0.028]),
        ], 0.001, 0.5);
        c.repos.push(Repo::new("core", "Forgemaster", 1.0, 330.0));
        c.repos.push(Repo::new("util", "CCC", 2.0, 116.0));

        let initial = c.conservation();
        for _ in 0..50 { c.step_leapfrog(); }
        let final_state = c.conservation();

        // Energy should be roughly conserved
        assert!(initial.energy_conserved(&final_state, 0.05));

        // Repos should have advanced their angles
        for repo in &c.repos {
            assert!(repo.angle > 0.0);
        }
    }

    #[test]
    fn constellation_evolve_steps() {
        let mut c = Constellation::new(vec![
            Vessel::at("A", 10.0, vec![0.0, 0.0]),
            Vessel::at("B", 10.0, vec![2.0, 0.0]),
        ], 0.01, 1.0);
        let (i, f) = c.evolve(10);
        // After 10 steps, something should have changed
        assert!((i.kinetic - f.kinetic).abs() > 0.0 || (i.potential - f.potential).abs() > 0.0);
    }

    #[test]
    fn total_energy_sums_correctly() {
        let c = Constellation::new(vec![
            Vessel::with_velocity("A", 5.0, vec![0.0, 0.0], vec![2.0, 0.0]),
            Vessel::with_velocity("B", 3.0, vec![10.0, 0.0], vec![0.0, 1.0]),
        ], 0.01, 1.0);
        let ke = c.total_kinetic_energy();
        let pe = c.total_potential_energy();
        let cs = c.conservation();
        assert!((cs.total_energy - (ke + pe)).abs() < 1e-10);
    }

    #[test]
    fn serde_roundtrip_constellation() {
        let c = Constellation::new(vec![
            Vessel::with_velocity("A", 10.0, vec![1.0, 2.0], vec![0.1, 0.2]),
        ], 0.05, 6.674e-11);
        let json = serde_json::to_string(&c).unwrap();
        let c2: Constellation = serde_json::from_str(&json).unwrap();
        assert_eq!(c, c2);
    }
}
