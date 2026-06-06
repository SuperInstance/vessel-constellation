# 🌌 vessel-constellation

**An N-body gravitational simulation of a fleet of software vessels and their orbiting repositories.**

```
          ★ Oracle (43)
         ╱
        ╱    ─ ─ ─ JetsonClaw1 (76)
       ╱   ╱
      ╱  ╱     ╭─────────────╮
     ╱ ╱    ┌──│─────────────│──┐
    ╱╱      │  │  Forgemaster │  │   Each dot = a repo
   ★        │  │   (330)      │  │   orbiting its vessel
    ╲╲      │  │              │  │
     ╲ ╲    └──│─────────────│──┘
      ╲  ╲     ╰─────────────╯
       ╲   ╲        │
        ╲    ─ ─ ─ ★ CCC (116)
         ╲
          ★
```

## The Metaphor

Imagine your software organization as a **solar system**:

- **Vessels** are stars — massive gravitational bodies whose weight comes from the number of repositories they contain
- **Repositories** are planets — orbiting their parent vessel, with core repos close and fast, peripheral repos far and slow
- **Dependencies** are gravity — pulling vessels toward each other with force proportional to the product of their masses and inversely proportional to the square of the distance between them
- **Conservation laws** govern the fleet: total energy and angular momentum are preserved as the constellation evolves

### The Fleet

| Vessel        | Repos | Role                    |
|---------------|-------|-------------------------|
| Forgemaster   | 330   | The titan — anchors the system |
| CCC           | 116   | The pillar — stabilizes the middle |
| JetsonClaw1   | 76    | The operative — bridges the gap |
| Oracle        | 43    | The sentinel — scouts the frontier |

```
  Mass → Gravitational Pull
  ━━━━━━━━━━━━━━━━━━━━━━━━
  Forgemaster (330) ████████████████████████████  ← dominates
  CCC          (116) ██████████
  JetsonClaw1   (76) ███████
  Oracle        (43) ████
```

## Architecture

```
┌─────────────────────────────────────────────────────┐
│                   Constellation                      │
│  ┌─────────┐  ┌─────────┐  ┌─────────┐  ┌────────┐ │
│  │Forgemstr│  │  CCC    │  │JetsonClw│  │ Oracle │ │
│  │ m=330   │  │ m=116   │  │ m=76    │  │ m=43   │ │
│  │    ●    │  │   ●     │  │  ●      │  │ ●      │ │
│  │  ● ● ● │  │ ● ● ●  │  │ ● ●    │  │ ●      │ │
│  │ ● ● ●● │  │● ● ●●  │  │● ● ●   │  │● ●     │ │
│  └─────────┘  └─────────┘  └─────────┘  └────────┘ │
│       ↕            ↕            ↕            ↕      │
│  ═════════ Gravitational Field (F = Gm₁m₂/r²) ═══  │
│       ↕            ↕            ↕            ↕      │
│              Conservation Laws                       │
│         E = KE + PE = const                         │
│         L = Σ m(r × v) = const                      │
└─────────────────────────────────────────────────────┘
```

## Modules

| Module            | Purpose                                              |
|-------------------|------------------------------------------------------|
| `vessel`          | Vessel as gravitational body: mass, position, velocity |
| `orbit`           | Repo orbits: Keplerian dynamics, T² ∝ r³            |
| `gravity`         | N-body forces, potential energy, Lagrange detection  |
| `conservation`    | Track total energy E and angular momentum L          |
| `perturbation`    | Events: repo added/removed, dependency shifts        |
| `constellation`   | Full state + leapfrog integration                    |

## Quick Start

```rust
use vessel_constellation::*;

// Create the fleet
let mut constellation = constellation::initial_fleet();

// Add orbiting repos
constellation.repos.push(Repo::new("core-engine", "Forgemaster", 1.0, 330.0));
constellation.repos.push(Repo::new("ml-pipeline", "JetsonClaw1", 3.0, 76.0));

// Evolve the system through 1000 time steps
let (initial, final_state) = constellation.evolve(1000);

// Verify conservation laws
assert!(initial.energy_conserved(&final_state, 0.01));
assert!(initial.angular_momentum_conserved(&final_state, 0.01));
```

## Gravitational Dynamics

### Force Between Vessels

```
  Forgemaster                    CCC
       ★ ←───────────────────── ★
              F = G·m₁·m₂/r²
         = 1.0 × 330 × 116 / 100
              = 382.8
```

### Orbital Mechanics

Repos orbit their vessel following Kepler's third law:

```
  ω = √(μ / r³)     angular velocity
  T = 2π / ω         orbital period
  T² = 4π²r³/μ       Kepler's 3rd law

  Core repo (r=1):     T = 0.35   ← fast, tight orbit
  Standard repo (r=3): T = 1.82
  Peripheral (r=10):   T = 11.07  ← slow, wide orbit
```

```
        ╭───────────────────────────╮
        │         ○ ← r=10          │  Peripheral repos
        │      ○    ●    ○          │
        │    ○   ○ ★ ○   ○         │  Standard repos
        │      ○  ● ●  ○           │
        │        ○← r=1            │  Core repos
        ╰───────────────────────────╯
                    ▲
              Vessel Center
```

### Conservation Laws

The leapfrog integrator is **symplectic** — it preserves the geometric structure of Hamiltonian mechanics:

```
  Total Energy:     E = ½Σmᵢvᵢ² - GΣmᵢmⱼ/rᵢⱼ
  Angular Momentum: L = Σmᵢ(rᵢ × vᵢ)

  ┌──────────────────────────────────────┐
  │  Step │     Energy  │  Ang. Momentum │
  ├───────┼─────────────┼────────────────┤
  │    0  │ -247.300000 │    50.000000   │
  │   50  │ -247.299981 │    50.000000   │
  │  100  │ -247.299963 │    50.000000   │  ← conserved!
  │  500  │ -247.299801 │    50.000000   │
  │ 1000  │ -247.299512 │    50.000000   │
  └───────┴─────────────┴────────────────┘
```

### Leapfrog vs Euler

```
  Energy drift over 200 steps:
  ┌────────────────────────────────────────┐
  │ Leapfrog: ▏         (drift < 0.001%)  │
  │ Euler:    ████████████ (drift ~5-15%)  │
  └────────────────────────────────────────┘
```

## Perturbation Events

```
  ┌──────────────┐     ┌──────────────┐
  │  Stable       │     │  Perturbed    │
  │  Constellation│────→│  Constellation│
  │              │     │              │
  │  FM  CCC  JC │     │  FM  CCC  JC │
  │  330 116  76 │     │ 331 116  76 │  ← +1 repo to FM
  │    O    O    │     │    O    O    │
  └──────────────┘     └──────────────┘
         │                     │
         │  ΔE = +0.04        │  Gravitational field
         │  ΔL = +0.12        │  adjusts to new mass
         └─────────────────────┘
```

```rust
use vessel_constellation::perturbation::Perturbation;

// Add a new repo to Forgemaster
let event = Perturbation::RepoAdded {
    vessel_name: "Forgemaster".into(),
    repo_name: "new-library".into(),
    mass_delta: 1.0,
};

// Apply and measure impact
let delta = event.conservation_delta(&vessels, &field);
println!("Energy change: {}", delta.total_energy);
```

## Lagrange Points

Detect when three vessels form a stable equilateral configuration:

```
        ★ JetsonClaw1
       ╱ ╲
      ╱   ╲
     ╱  L₄  ╲
    ╱       ╲
   ★─────────★
  Forgemaster  CCC
        L₅
```

```rust
use vessel_constellation::gravity::is_lagrange_triangle;

// Equilateral triangle → stable Lagrange configuration
let a = Vessel::at("A", 10.0, vec![0.0, 0.0]);
let b = Vessel::at("B", 10.0, vec![2.0, 0.0]);
let c = Vessel::at("C", 10.0, vec![1.0, 3.0_f64.sqrt()]);
assert!(is_lagrange_triangle(&a, &b, &c));
```

## API Reference

### Core Types

```rust
struct Vessel {
    name: String,
    mass: f64,              // = repo count
    position: Vec<f64>,     // dependency-space coordinates
    velocity: Vec<f64>,     // growth rate vector
}

struct Repo {
    name: String,
    vessel: String,
    orbital_radius: f64,    // distance from vessel center
    angle: f64,             // current orbital angle
    angular_velocity: f64,  // from Kepler: ω = √(μ/r³)
}

struct Constellation {
    vessels: Vec<Vessel>,
    repos: Vec<Repo>,
    dt: f64,               // time step
    G: f64,                // gravitational constant
}
```

### Key Methods

```rust
// Vessel
vessel.distance_to(&other)      // Euclidean distance
vessel.kinetic_energy()         // ½mv²
vessel.angular_momentum()       // m(r × v)

// Repo
repo.period()                   // T = 2π/ω
repo.position()                 // (x, y) cartesian
repo.verify_kepler(mu)          // check T²μ = 4π²r³

// Constellation
constellation.step_leapfrog()   // symplectic integration
constellation.evolve(100)       // returns (initial, final) conservation states
constellation.conservation()    // current E, L, KE, PE
```

## Installation

```toml
[dependencies]
vessel-constellation = "0.1"
```

## License

MIT
