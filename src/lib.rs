//! # Vessel Constellation
//!
//! An N-body gravitational simulation where software vessels and their
//! orbiting repositories form a celestial system governed by Newtonian
//! mechanics, conservation laws, and Keplerian orbital dynamics.

#![allow(non_snake_case)]

pub mod vessel;
pub mod orbit;
pub mod gravity;
pub mod conservation;
pub mod perturbation;
pub mod constellation;

pub use vessel::Vessel;
pub use orbit::Repo;
pub use gravity::GravitationalField;
pub use conservation::ConservationState;
pub use constellation::Constellation;
