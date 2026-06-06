//! Conservation laws: total energy and angular momentum.
//!
//! L = Σ mᵢ(rᵢ × vᵢ) is conserved.
//! E = KE + PE is conserved.
//! Both should remain constant (within numerical tolerance) during
//! time evolution.

use serde::{Deserialize, Serialize};
use crate::vessel::Vessel;
use crate::gravity::GravitationalField;

/// Snapshot of conserved quantities.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ConservationState {
    pub total_energy: f64,
    pub angular_momentum: Vec<f64>,
    pub kinetic: f64,
    pub potential: f64,
}

impl ConservationState {
    /// Compute conservation state from vessels and gravitational field.
    pub fn compute(vessels: &[Vessel], field: &GravitationalField) -> Self {
        let kinetic: f64 = vessels.iter().map(|v| v.kinetic_energy()).sum();
        let potential = field.total_potential_energy();
        let total_energy = kinetic + potential;

        let angular_momentum = vessels.iter()
            .fold(vec![0.0, 0.0, 0.0], |acc, v| {
                let lam = v.angular_momentum();
                acc.iter().zip(lam.iter()).map(|(a, l)| a + l).collect()
            });

        Self { total_energy, angular_momentum, kinetic, potential }
    }

    /// Check if total energy is conserved relative to another state.
    pub fn energy_conserved(&self, other: &ConservationState, tolerance: f64) -> bool {
        let e_ref = self.total_energy.abs().max(1.0);
        (self.total_energy - other.total_energy).abs() / e_ref < tolerance
    }

    /// Check if angular momentum is conserved relative to another state.
    pub fn angular_momentum_conserved(&self, other: &ConservationState, tolerance: f64) -> bool {
        for (a, b) in self.angular_momentum.iter().zip(other.angular_momentum.iter()) {
            let ref_val = a.abs().max(1.0);
            if (a - b).abs() / ref_val > tolerance {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_two_body() -> (Vec<Vessel>, GravitationalField) {
        let G = 1.0;
        let vessels = vec![
            Vessel::with_velocity("A", 10.0, vec![0.0, 0.0], vec![0.0, 0.0]),
            Vessel::with_velocity("B", 10.0, vec![5.0, 0.0], vec![0.0, 2.0]),
        ];
        let field = GravitationalField::new(vessels.clone(), G);
        (vessels, field)
    }

    #[test]
    fn conservation_state_computation() {
        let (vessels, field) = make_two_body();
        let state = ConservationState::compute(&vessels, &field);
        assert!(state.kinetic > 0.0);
        assert!(state.potential < 0.0);
    }

    #[test]
    fn energy_components_sum() {
        let (vessels, field) = make_two_body();
        let state = ConservationState::compute(&vessels, &field);
        assert!((state.total_energy - (state.kinetic + state.potential)).abs() < 1e-10);
    }

    #[test]
    fn energy_conserved_check() {
        let s1 = ConservationState {
            total_energy: -100.0, angular_momentum: vec![0.0, 0.0, 50.0],
            kinetic: 20.0, potential: -120.0,
        };
        let s2 = ConservationState {
            total_energy: -100.001, angular_momentum: vec![0.0, 0.0, 50.001],
            kinetic: 20.001, potential: -120.002,
        };
        assert!(s1.energy_conserved(&s2, 1e-3));
    }

    #[test]
    fn angular_momentum_conserved_check() {
        let s1 = ConservationState {
            total_energy: -50.0, angular_momentum: vec![0.0, 0.0, 100.0],
            kinetic: 10.0, potential: -60.0,
        };
        let s2 = ConservationState {
            total_energy: -50.0, angular_momentum: vec![0.0, 0.0, 100.01],
            kinetic: 10.0, potential: -60.0,
        };
        assert!(s1.angular_momentum_conserved(&s2, 1e-3));
    }

    #[test]
    fn serde_roundtrip_conservation() {
        let s = ConservationState {
            total_energy: -42.5, angular_momentum: vec![0.0, 0.0, 7.3],
            kinetic: 10.0, potential: -52.5,
        };
        let json = serde_json::to_string(&s).unwrap();
        let s2: ConservationState = serde_json::from_str(&json).unwrap();
        assert_eq!(s, s2);
    }
}
