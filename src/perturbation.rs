//! Perturbation events that disturb the constellation.
//!
//! New repo added, repo deleted, dependency shift — each perturbation
//! propagates through the N-body system, altering masses, positions,
//! and velocities.

use serde::{Deserialize, Serialize};
use crate::vessel::Vessel;
use crate::orbit::Repo;
use crate::gravity::GravitationalField;
use crate::conservation::ConservationState;

/// A perturbation event affecting the constellation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum Perturbation {
    /// New repository added to a vessel. Increases mass.
    RepoAdded { vessel_name: String, repo_name: String, mass_delta: f64 },
    /// Repository removed from a vessel. Decreases mass.
    RepoRemoved { vessel_name: String, repo_name: String, mass_delta: f64 },
    /// Dependency shift changes a vessel's position in dependency-space.
    DependencyShift { vessel_name: String, delta_position: Vec<f64> },
    /// Velocity kick (e.g., sudden growth spurt).
    VelocityKick { vessel_name: String, delta_velocity: Vec<f64> },
}

impl Perturbation {
    /// Apply this perturbation to a mutable slice of vessels.
    /// Returns true if applied successfully.
    pub fn apply(&self, vessels: &mut [Vessel], repos: &mut Vec<Repo>) -> bool {
        match self {
            Perturbation::RepoAdded { vessel_name, repo_name, mass_delta } => {
                if let Some(v) = vessels.iter_mut().find(|v| v.name == *vessel_name) {
                    v.mass += mass_delta;
                    let mu = v.mass; // Use new mass as μ
                    repos.push(Repo::new(repo_name, vessel_name, 1.0, mu));
                    true
                } else { false }
            }
            Perturbation::RepoRemoved { vessel_name, repo_name, mass_delta } => {
                if let Some(v) = vessels.iter_mut().find(|v| v.name == *vessel_name) {
                    v.mass = (v.mass - mass_delta).max(0.0);
                    repos.retain(|r| r.name != *repo_name);
                    true
                } else { false }
            }
            Perturbation::DependencyShift { vessel_name, delta_position } => {
                if let Some(v) = vessels.iter_mut().find(|v| v.name == *vessel_name) {
                    for (p, d) in v.position.iter_mut().zip(delta_position.iter()) {
                        *p += d;
                    }
                    true
                } else { false }
            }
            Perturbation::VelocityKick { vessel_name, delta_velocity } => {
                if let Some(v) = vessels.iter_mut().find(|v| v.name == *vessel_name) {
                    for (vel, d) in v.velocity.iter_mut().zip(delta_velocity.iter()) {
                        *vel += d;
                    }
                    true
                } else { false }
            }
        }
    }

    /// Compute the change in conservation state due to this perturbation.
    pub fn conservation_delta(&self, vessels: &[Vessel], field: &GravitationalField) -> ConservationState {
        let before = ConservationState::compute(vessels, field);
        let mut vessels_copy = vessels.to_vec();
        let mut repos_copy: Vec<Repo> = Vec::new();
        self.apply(&mut vessels_copy, &mut repos_copy);
        let new_field = GravitationalField::new(vessels_copy.clone(), field.G);
        let after = ConservationState::compute(&vessels_copy, &new_field);

        ConservationState {
            total_energy: after.total_energy - before.total_energy,
            angular_momentum: after.angular_momentum.iter()
                .zip(before.angular_momentum.iter())
                .map(|(a, b)| a - b)
                .collect(),
            kinetic: after.kinetic - before.kinetic,
            potential: after.potential - before.potential,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repo_added_increases_mass() {
        let mut vessels = vec![Vessel::new("Forgemaster", 330.0, 2)];
        let mut repos = Vec::new();
        let p = Perturbation::RepoAdded { vessel_name: "Forgemaster".into(), repo_name: "new-lib".into(), mass_delta: 1.0 };
        p.apply(&mut vessels, &mut repos);
        assert_eq!(vessels[0].mass, 331.0);
        assert_eq!(repos.len(), 1);
    }

    #[test]
    fn repo_removed_decreases_mass() {
        let mut vessels = vec![Vessel::new("CCC", 116.0, 2)];
        let mut repos = vec![Repo::new("old-lib", "CCC", 1.0, 116.0)];
        let p = Perturbation::RepoRemoved { vessel_name: "CCC".into(), repo_name: "old-lib".into(), mass_delta: 1.0 };
        p.apply(&mut vessels, &mut repos);
        assert_eq!(vessels[0].mass, 115.0);
        assert!(repos.is_empty());
    }

    #[test]
    fn dependency_shift_moves_position() {
        let mut vessels = vec![Vessel::at("JetsonClaw1", 76.0, vec![1.0, 2.0])];
        let mut repos = Vec::new();
        let p = Perturbation::DependencyShift { vessel_name: "JetsonClaw1".into(), delta_position: vec![0.5, -0.3] };
        p.apply(&mut vessels, &mut repos);
        assert!((vessels[0].position[0] - 1.5).abs() < 1e-10);
        assert!((vessels[0].position[1] - 1.7).abs() < 1e-10);
    }

    #[test]
    fn velocity_kick_changes_velocity() {
        let mut vessels = vec![Vessel::with_velocity("Oracle", 43.0, vec![0.0, 0.0], vec![0.1, 0.0])];
        let mut repos = Vec::new();
        let p = Perturbation::VelocityKick { vessel_name: "Oracle".into(), delta_velocity: vec![0.0, 0.5] };
        p.apply(&mut vessels, &mut repos);
        assert!((vessels[0].velocity[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn perturbation_on_missing_vessel_fails() {
        let mut vessels = vec![Vessel::new("A", 10.0, 2)];
        let mut repos = Vec::new();
        let p = Perturbation::RepoAdded { vessel_name: "Z".into(), repo_name: "x".into(), mass_delta: 1.0 };
        assert!(!p.apply(&mut vessels, &mut repos));
    }

    #[test]
    fn conservation_delta_tracks_energy_change() {
        let vessels = vec![
            Vessel::at("A", 10.0, vec![0.0, 0.0]),
            Vessel::at("B", 10.0, vec![5.0, 0.0]),
        ];
        let field = GravitationalField::new(vessels.clone(), 1.0);
        let p = Perturbation::VelocityKick { vessel_name: "A".into(), delta_velocity: vec![1.0, 0.0] };
        let delta = p.conservation_delta(&vessels, &field);
        // Adding velocity should increase kinetic energy
        assert!(delta.kinetic > 0.0);
    }

    #[test]
    fn serde_roundtrip_perturbation() {
        let p = Perturbation::DependencyShift {
            vessel_name: "Forgemaster".into(),
            delta_position: vec![0.1, -0.2, 0.3],
        };
        let json = serde_json::to_string(&p).unwrap();
        let p2: Perturbation = serde_json::from_str(&json).unwrap();
        assert_eq!(p, p2);
    }
}
