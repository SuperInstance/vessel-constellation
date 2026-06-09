//! # vessel-constellation Tutorial
//!
//! Learn how software vessels and repositories form an N-body gravitational system.
//!
//! **Mathematical insight:** Four vessels (Forgemaster, CCC, JetsonClaw1, Oracle)
//! orbit each other under Newtonian gravity F = G·m₁·m₂/r². Repositories orbit
//! their vessel following Kepler's third law: T² ∝ r³. The system is integrated
//! using symplectic leapfrog (Störmer-Verlet) for superior energy conservation.
//!
//! Conservation laws (energy E, angular momentum L) are the invariant backbone.
//!
//! Run: `cargo run --example tutorial`

use vessel_constellation::{
    Vessel, Repo, GravitationalField, ConservationState,
    Constellation,
    perturbation::Perturbation,
};

fn main() {
    println!("════════════════════════════════════════════════════════");
    println!("  vessel-constellation: N-Body Orbital Dynamics        ");
    println!("════════════════════════════════════════════════════════\n");

    lesson_1_vessels_as_bodies();
    lesson_2_gravitational_forces();
    lesson_3_kepler_orbits();
    lesson_4_conservation_laws();
    lesson_5_leapfrog_integration();
    lesson_6_perturbations();
    lesson_7_full_constellation();
}

/// Lesson 1: Vessels as Gravitational Bodies.
///
/// Each vessel is a point mass in dependency-space:
/// - Mass = number of repositories (proxy for gravitational pull)
/// - Position = location in coupling-category space
/// - Velocity = growth rate vector
///
/// Forgemaster (330 repos) dominates. Oracle (43 repos) is a minor body.
fn lesson_1_vessels_as_bodies() {
    println!("━━━ Lesson 1: Vessels as Gravitational Bodies ━━━\n");

    let forgemaster = Vessel::at("Forgemaster", 330.0, vec![0.0, 0.0]);
    let ccc = Vessel::at("CCC", 116.0, vec![10.0, 0.0]);
    let jetson = Vessel::at("JetsonClaw1", 76.0, vec![5.0, 8.0]);
    let oracle = Vessel::at("Oracle", 43.0, vec![-3.0, 5.0]);

    let vessels = [&forgemaster, &ccc, &jetson, &oracle];
    println!("  Fleet composition:");
    for v in &vessels {
        println!("    {:>12}: mass={:>6.0}, position=({:>5.1}, {:>5.1})",
                 v.name, v.mass, v.position[0], v.position[1]);
    }

    // Distance matrix
    println!("\n  Distance matrix:");
    print!("              ");
    for v in &vessels {
        print!("{:>12}", v.name.chars().take(8).collect::<String>());
    }
    println!();
    for va in &vessels {
        print!("  {:>12} ", va.name.chars().take(10).collect::<String>());
        for vb in &vessels {
            let d = va.distance_to(vb);
            print!("{:>12.2}", d);
        }
        println!();
    }

    // Center of mass
    let com = forgemaster.center_of_mass(&ccc);
    println!("\n  Center of mass (Forgemaster + CCC): ({:.2}, {:.2})", com[0], com[1]);
    println!();
}

/// Lesson 2: Gravitational Forces — F = G·m₁·m₂/r².
///
/// Newton's law of gravitation: force is attractive, proportional to
/// product of masses, inversely proportional to distance squared.
/// Heavier vessels pull harder. Closer vessels pull harder.
fn lesson_2_gravitational_forces() {
    println!("━━━ Lesson 2: Gravitational Forces ━━━\n");

    let fm = Vessel::at("Forgemaster", 330.0, vec![0.0, 0.0]);
    let ccc = Vessel::at("CCC", 116.0, vec![10.0, 0.0]);
    let oracle = Vessel::at("Oracle", 43.0, vec![-3.0, 5.0]);

    let field = GravitationalField::new(vec![fm.clone(), ccc.clone(), oracle.clone()], 1.0);

    // Force between each pair
    println!("  Pairwise forces (G=1.0):");
    let pairs = [(&fm, &ccc), (&fm, &oracle), (&ccc, &oracle)];
    for (a, b) in &pairs {
        let force = field.force_between(a, b);
        let magnitude = force.iter().map(|f| f * f).sum::<f64>().sqrt();
        println!("    {} → {}: F=({:.3}, {:.3}) |F|={:.3}",
                 a.name, b.name, force[0], force[1], magnitude);
    }

    // Net force on Oracle (being pulled by both big vessels)
    let net_oracle = field.net_force(&oracle);
    let net_mag = net_oracle.iter().map(|f| f * f).sum::<f64>().sqrt();
    println!("\n  Net force on Oracle: ({:.3}, {:.3}), |F|={:.3}",
             net_oracle[0], net_oracle[1], net_mag);

    // Potential energy
    let pe = field.total_potential_energy();
    println!("\n  Total potential energy: {:.3} (negative = bound system)", pe);

    // Gravitational potential at a point
    let phi = field.potential_at(&[5.0, 5.0]);
    println!("  Gravitational potential at (5,5): {:.3}", phi);

    // Accelerations (F/m — lighter bodies accelerate faster)
    let accs = field.accelerations();
    println!("\n  Accelerations:");
    for (v, a) in field.sources.iter().zip(accs.iter()) {
        let a_mag = a.iter().map(|x| x * x).sum::<f64>().sqrt();
        println!("    {:>12}: a=({:.4}, {:.4}) |a|={:.4}", v.name, a[0], a[1], a_mag);
    }
    println!();
}

/// Lesson 3: Kepler Orbits — T² ∝ r³.
///
/// Repositories orbit their vessel following Kepler-like dynamics:
///   ω = √(μ/r³)  where μ = G × vessel_mass
///   T = 2π/ω = 2π√(r³/μ)
///
/// Core repos (small r) orbit fast. Peripheral repos (large r) orbit slow.
/// This is the same law that governs planets around a star.
fn lesson_3_kepler_orbits() {
    println!("━━━ Lesson 3: Kepler Orbits — T² ∝ r³ ━━━\n");

    let mu = 330.0; // G × Forgemaster mass

    // Create repos at different orbital radii
    let core = Repo::new("core-engine", "Forgemaster", 1.0, mu);
    let mid = Repo::new("web-framework", "Forgemaster", 3.0, mu);
    let outer = Repo::new("legacy-tool", "Forgemaster", 8.0, mu);

    let repos = [&core, &mid, &outer];
    println!("  Forgemaster orbital mechanics (μ={}):\n", mu);

    for r in &repos {
        let (x, y) = r.position();
        println!("    {:>15}: r={:.1}, ω={:.4}, T={:.2}, v={:.2}, E={:.2}",
                 r.name, r.orbital_radius, r.angular_velocity, r.period(),
                 r.orbital_speed(), r.orbital_energy(mu));
        println!("      initial position: ({:.2}, {:.2})", x, y);
    }

    // Verify Kepler's third law
    println!("\n  Kepler's third law verification:");
    for r in &repos {
        let valid = r.verify_kepler(mu);
        println!("    {}: T²·μ = {:.2}, 4π²r³ = {:.2}, valid={}",
                 r.name, r.period().powi(2) * mu,
                 4.0 * std::f64::consts::PI.powi(2) * r.orbital_radius.powi(3),
                 valid);
    }

    // Advance orbits and show positions
    println!("\n  Core repo positions over time:");
    let mut sim_core = Repo::new("core-engine", "Forgemaster", 1.0, mu);
    for step in 0..8 {
        let (x, y) = sim_core.position();
        println!("    t={:.2}: ({:>6.3}, {:>6.3})  angle={:.2} rad",
                 step as f64 * 0.5, x, y, sim_core.angle);
        sim_core.step(0.5);
    }
    println!();
}

/// Lesson 4: Conservation Laws — E and L are invariant.
///
/// In a closed gravitational system:
/// - Total energy E = KE + PE is conserved
/// - Total angular momentum L = Σ mᵢ(rᵢ × vᵢ) is conserved
///
/// These are the deep symmetries: energy conservation from time-translation
/// invariance, angular momentum from rotational invariance (Noether's theorem).
fn lesson_4_conservation_laws() {
    println!("━━━ Lesson 4: Conservation Laws ━━━\n");

    // Two-body system for clarity
    let vessels = vec![
        Vessel::with_velocity("A", 10.0, vec![0.0, 0.0], vec![0.0, 2.0]),
        Vessel::with_velocity("B", 10.0, vec![5.0, 0.0], vec![0.0, -2.0]),
    ];
    let field = GravitationalField::new(vessels.clone(), 1.0);

    let state = ConservationState::compute(&vessels, &field);
    println!("  Initial conservation state:");
    println!("    Total energy:   {:.4}", state.total_energy);
    println!("    Kinetic energy: {:.4}", state.kinetic);
    println!("    Potential energy: {:.4}", state.potential);
    println!("    Angular momentum: {:?}", state.angular_momentum);

    // Verify E = KE + PE
    let e_check = state.kinetic + state.potential;
    println!("\n  E = KE + PE check: {:.4} = {:.4} + {:.4} = {:.4} ✓",
             state.total_energy, state.kinetic, state.potential, e_check);

    // Compare two states
    let state2 = ConservationState {
        total_energy: state.total_energy + 0.001,
        angular_momentum: state.angular_momentum.clone(),
        kinetic: state.kinetic + 0.001,
        potential: state.potential,
    };
    println!("\n  Energy conserved (tolerance 0.01)? {}", state.energy_conserved(&state2, 0.01));
    println!("  Angular momentum conserved? {}", state.angular_momentum_conserved(&state2, 0.01));
    println!();
}

/// Lesson 5: Leapfrog Integration — Symplectic beats Euler.
///
/// The leapfrog (Störmer-Verlet) integrator is symplectic: it preserves
/// the geometric structure of Hamiltonian mechanics. Result: energy
/// oscillates around the true value instead of drifting.
///
/// Comparison:
///   Euler: E drifts linearly → simulation becomes unphysical
///   Leapfrog: E oscillates → long-term stability
fn lesson_5_leapfrog_integration() {
    println!("━━━ Lesson 5: Leapfrog vs Euler Integration ━━━\n");

    let make_system = || Constellation::new(vec![
        Vessel::with_velocity("A", 10.0, vec![0.0, 0.0], vec![0.0, 0.3]),
        Vessel::with_velocity("B", 10.0, vec![4.0, 0.0], vec![0.0, -0.3]),
    ], 0.005, 1.0);

    // Leapfrog
    let mut lf = make_system();
    let lf_initial = lf.conservation();
    for _ in 0..500 { lf.step_leapfrog(); }
    let lf_final = lf.conservation();
    let lf_drift = (lf_final.total_energy - lf_initial.total_energy).abs();

    // Euler
    let mut eu = make_system();
    let eu_initial = eu.conservation();
    for _ in 0..500 { eu.step_euler(); }
    let eu_final = eu.conservation();
    let eu_drift = (eu_final.total_energy - eu_initial.total_energy).abs();

    println!("  500 steps, dt=0.005:");
    println!("    Leapfrog energy drift: {:.6}", lf_drift);
    println!("    Euler energy drift:    {:.6}", eu_drift);
    println!("    Leapfrog is {:.1}× better", eu_drift / lf_drift.max(1e-15));

    // Show conservation over time with leapfrog
    println!("\n  Leapfrog energy conservation over time:");
    let mut c = make_system();
    let initial_e = c.conservation().total_energy;
    for step in [0, 100, 200, 300, 400, 500] {
        let state = c.conservation();
        let drift = (state.total_energy - initial_e).abs();
        println!("    step {:>4}: E={:>10.4}, drift={:.6}", step, state.total_energy, drift);
        for _ in 0..100 { c.step_leapfrog(); }
    }

    assert!(lf_drift < eu_drift, "Leapfrog should conserve energy better");
    println!();
}

/// Lesson 6: Perturbations — Disturbing the constellation.
///
/// Real systems aren't static. New repos are added, old ones removed,
/// dependencies shift. Each perturbation changes the gravitational field
/// and propagates through the constellation.
///
/// Conservation laws are temporarily violated by perturbations —
/// the system must re-equilibrate.
fn lesson_6_perturbations() {
    println!("━━━ Lesson 6: Perturbation Events ━━━\n");

    let mut fleet = Constellation::new(vec![
        Vessel::at("Forgemaster", 330.0, vec![0.0, 0.0]),
        Vessel::at("CCC", 116.0, vec![10.0, 0.0]),
    ], 0.001, 1.0);

    let initial = fleet.conservation();
    println!("  Initial state: E={:.4}", initial.total_energy);

    // Perturbation 1: Add a repo to CCC
    let p1 = Perturbation::RepoAdded {
        vessel_name: "CCC".into(),
        repo_name: "new-lib".into(),
        mass_delta: 5.0,
    };
    let field = fleet.field();
    let delta1 = p1.conservation_delta(&fleet.vessels, &field);
    p1.apply(&mut fleet.vessels, &mut fleet.repos);
    println!("\n  After adding 'new-lib' to CCC (+5 mass):");
    println!("    CCC mass: {}", fleet.vessels[1].mass);
    println!("    Energy delta: {:.4} (KE={:.4}, PE={:.4})",
             delta1.total_energy, delta1.kinetic, delta1.potential);

    // Perturbation 2: Dependency shift
    let p2 = Perturbation::DependencyShift {
        vessel_name: "Forgemaster".into(),
        delta_position: vec![0.5, -0.3],
    };
    let field2 = fleet.field();
    let delta2 = p2.conservation_delta(&fleet.vessels, &field2);
    p2.apply(&mut fleet.vessels, &mut fleet.repos);
    println!("\n  After dependency shift on Forgemaster:");
    println!("    Forgemaster position: ({:.1}, {:.1})",
             fleet.vessels[0].position[0], fleet.vessels[0].position[1]);
    println!("    Energy delta: {:.4}", delta2.total_energy);

    // Perturbation 3: Velocity kick (growth spurt)
    let p3 = Perturbation::VelocityKick {
        vessel_name: "CCC".into(),
        delta_velocity: vec![0.0, 0.1],
    };
    p3.apply(&mut fleet.vessels, &mut fleet.repos);
    println!("\n  After velocity kick on CCC:");
    println!("    CCC velocity: ({:.3}, {:.3})",
             fleet.vessels[1].velocity[0], fleet.vessels[1].velocity[1]);

    // Perturbation 4: Remove a repo
    let p4 = Perturbation::RepoRemoved {
        vessel_name: "CCC".into(),
        repo_name: "new-lib".into(),
        mass_delta: 5.0,
    };
    p4.apply(&mut fleet.vessels, &mut fleet.repos);
    println!("\n  After removing 'new-lib' from CCC:");
    println!("    CCC mass: {}, repos remaining: {}", fleet.vessels[1].mass, fleet.repos.len());
    println!();
}

/// Lesson 7: Full Constellation Simulation.
///
/// Put it all together: four vessels with orbiting repos, evolving under
/// gravity with leapfrog integration. Watch conservation laws hold and
/// the system settle into orbital dynamics.
fn lesson_7_full_constellation() {
    println!("━━━ Lesson 7: Full Constellation Evolution ━━━\n");

    let mut constellation = Constellation::new(vec![
        Vessel::with_velocity("Forgemaster", 330.0, vec![0.0, 0.0], vec![0.0, 0.01]),
        Vessel::with_velocity("CCC", 116.0, vec![8.0, 0.0], vec![0.0, -0.028]),
        Vessel::with_velocity("JetsonClaw1", 76.0, vec![5.0, 8.0], vec![-0.01, 0.0]),
        Vessel::with_velocity("Oracle", 43.0, vec![-3.0, 5.0], vec![0.015, 0.0]),
    ], 0.001, 0.5);

    // Add orbiting repos
    constellation.repos.push(Repo::new("core", "Forgemaster", 0.5, 0.5 * 330.0));
    constellation.repos.push(Repo::new("stdlib", "Forgemaster", 1.5, 0.5 * 330.0));
    constellation.repos.push(Repo::new("toolkit", "CCC", 1.0, 0.5 * 116.0));

    let initial = constellation.conservation();
    println!("  Initial fleet:");
    for v in &constellation.vessels {
        println!("    {:>12}: mass={:>6.0}, pos=({:>6.2}, {:>6.2}), vel=({:>7.4}, {:>7.4})",
                 v.name, v.mass, v.position[0], v.position[1],
                 v.velocity[0], v.velocity[1]);
    }
    println!("  Initial energy: {:.4}", initial.total_energy);

    // Evolve for 500 steps
    println!("\n  Evolving 500 steps (dt=0.001)...");
    let (init_state, final_state) = constellation.evolve(500);

    println!("\n  Final fleet:");
    for v in &constellation.vessels {
        println!("    {:>12}: pos=({:>6.2}, {:>6.2}), vel=({:>7.4}, {:>7.4})",
                 v.name, v.position[0], v.position[1],
                 v.velocity[0], v.velocity[1]);
    }

    // Conservation check
    let energy_drift = (final_state.total_energy - init_state.total_energy).abs();
    let energy_ref = init_state.total_energy.abs().max(1.0);
    println!("\n  Conservation check:");
    println!("    Initial E: {:.6}", init_state.total_energy);
    println!("    Final E:   {:.6}", final_state.total_energy);
    println!("    Relative drift: {:.6} ({:.3}%)",
             energy_drift / energy_ref,
             energy_drift / energy_ref * 100.0);
    println!("    Energy conserved (1% tol): {}",
             init_state.energy_conserved(&final_state, 0.01));

    // Repo positions
    println!("\n  Orbiting repos after evolution:");
    for repo in &constellation.repos {
        let (x, y) = repo.position();
        println!("    {:>10} → {:>12}: pos=({:>6.3}, {:>6.3}), angle={:.2} rad",
                 repo.name, repo.vessel, x, y, repo.angle);
    }

    // Circular orbit detection
    use vessel_constellation::gravity::is_circular_orbit;
    let v0 = &constellation.vessels[0];
    let v1 = &constellation.vessels[1];
    let circular = is_circular_orbit(v0, v1, constellation.G);
    println!("\n  Forgemaster-CCC circular orbit? {}", circular);

    // Lagrange triangle check
    use vessel_constellation::gravity::is_lagrange_triangle;
    let lagrange = is_lagrange_triangle(
        &constellation.vessels[0],
        &constellation.vessels[1],
        &constellation.vessels[2],
    );
    println!("  FM-CCC-JC Lagrange triangle? {}", lagrange);

    println!("\n  ✦ Key insight: Symplectic leapfrog integration preserves the");
    println!("    Hamiltonian structure — energy oscillates rather than drifts.");
    println!("    This is why orbital mechanics simulations use leapfrog, not Euler.\n");
}
