use wasm_bindgen::prelude::*;
use std::ops::{Index, IndexMut};

/// Sets the global allocator to use 'wee_alloc'
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
extern {
    // External JavaScript prototypes
}

/// Setup some debugging utilities and any other onetime initialization
#[wasm_bindgen]
pub fn init() {

    #[cfg(feature = "console_error_panic_hook")]
    std::panic::set_hook(Box::new(console_error_panic_hook::hook));

    #[cfg(feature = "wasm-logger")]
    wasm_logger::init(wasm_logger::Config::default());
}

/// Universal Constants
const PERMITTIVITY: f32 = 8.8541878e-12; // Vacuum Permittivity [C/(V*m)]
const ELECTRON_CHARGE: f32 = 1.602176565e-19; // Electron Charge [C]
const AMU: f32 = 1.660538921e-27; // Atomic Mass Unit [kg]
const ELECTRON_MASS: f32 = 9.10938215e-31; // Electron Mass [kg]
const BOLTZMANN: f32 = 1.380648e-23; // Boltzmann Constant [J/K]

/*
 * init:
 * - initialize field
 * - solve initial potential
 * - compute the electric field
 * - load particles
 * 
 * loop:
 * - compute charge density
 * - compute potential
 * - compute electric field
 *
 * - update particle velocity
 * - update particle position
 * 
 * - perform collisions
 * - add / remove particles
 * - diagnostics
 *
 *
 * - output final results
 */ 

/// Simple Simulation Setup
#[wasm_bindgen]
pub fn sim() {

    let size = [10, 10, 10];

    let mut domain = Domain::new(size, [0.2, 0.2, 0.2]);

    log::info!("Node Volume: {}", domain.node_volume());

    let mut phi = Field::new(size);

    log::info!("Element: {}", phi[0][0][0]);

    let b = &mut phi[0][0][0];
    *b = 3.14;

    log::info!("Element: {}", phi[0][0][0]);

    let mut solver = PotentialSolver::new(&mut domain, 50, 0.001); 
    log::info!("Solver: {}", solver.solve());
}

/// Computes the Debye length in meters
///
/// * `temp` - Species temperature in K
/// * `density` - Species density in $$ \frac{#}{m^3} $$
fn debye_length(temp: f32, density: f32) -> f32 {

    // $$ \lambda_D = \sqrt{\frac{\epsilon_0 k_b T_e}{n_e {q_e}^2}} $$
    return ((PERMITTIVITY * BOLTZMANN * temp)/(density * ELECTRON_CHARGE * ELECTRON_CHARGE)).sqrt();
}

/// Represents the computational domain of the simulation as a Cartesian mesh
struct Domain {
    nodes: [usize; 3], // The number of nodes in each axis
    spacing: [f32; 3], // Space between nodes in x, y, and z
    phi: Field, // Electric Potential
    rho: Field, // Charge Density
    electric_field: Field, // Electric Field Components (TODO: This needs a field that stores a 3Vec)
}

impl Domain {

    /// Create a new domain
    pub fn new(nodes: [usize; 3], size: [f32; 3]) -> Domain {
        Domain {
            nodes: nodes,
            spacing: [
                size[0] / (nodes[0] - 1) as f32,
                size[1] / (nodes[1] - 1) as f32,
                size[2] / (nodes[2] - 1) as f32
            ],
            phi: Field::new(nodes),
            rho: Field::new(nodes),
            electric_field: Field::new(nodes)
        }
    }

    /// Compute the volumn of a cell
    pub fn node_volume(&self) -> f32 {
        self.spacing[0] * self.spacing[1] * self.spacing[2]
    }
}

/// 3D field of floats
/// TODO: Use generics to allow it to store arbitrary types
struct Field {
    data: Vec<Vec<Vec<f32>>>
}

impl Field {

    /// Create a new 3D field with dimensions x, y, z
    pub fn new(size: [usize; 3]) -> Field {
        Field {
            data: vec![vec![vec![0.0; size[2]]; size[1]]; size[0]]
        }
    }

    /// Create a new 3D field with equal dimensions
    pub fn new_cube(size: usize) -> Field {
        Field::new([size, size, size])
    }
}

/// Support indexing into a Field as if it were a Vec
impl Index<usize> for Field {
    type Output = Vec<Vec<f32>>;

    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}

/// Support mutable indexing into a Field as if it were a Vec
impl IndexMut<usize> for Field {
    fn index_mut(&mut self, i: usize) -> &mut Vec<Vec<f32>> {
        self.data.index_mut(i)
    }
}

struct Species {
    mass: f32,
    charge: f32,
    // Particle positions and velocities in the domain
    // TODO: Figure out what that looks like
}

struct PotentialSolver<'a> {
    domain: &'a mut Domain,
    iterations: usize,
    tolerance: f32
}

impl PotentialSolver<'_> {

    /// Create a new solver
    pub fn new(domain: &mut Domain, iterations: usize, tolerance: f32) -> PotentialSolver {
        PotentialSolver {
            domain,
            iterations,
            tolerance
        }
    }

    /// Iteratively computes the electrical potential at each node
    pub fn solve(&mut self) -> bool {

        // Get local references for ease of use
        let phi = &mut self.domain.phi;
        let rho = &self.domain.rho;

        // Precompute $$\frac{1}{\text{dx}^2}$$
        let s = self.domain.spacing;
        let dx2 = 1.0 / (s[0] * s[0]);
        let dy2 = 1.0 / (s[1] * s[1]);
        let dz2 = 1.0 / (s[2] * s[2]);

        // Solve the potential
        for iteration in 0..self.iterations {

            // Iterate over every internal element, since we are assuming 
            // the dirichlet boundary condition we don't process the edge elements
            for i in 1..self.domain.nodes[0]-1 {
                for j in 1..self.domain.nodes[1]-1 {
                    for k in 1..self.domain.nodes[2]-1 {

                        // Standard internal open node
                        let phi_new = ((rho[i][j][k]/PERMITTIVITY) + 
                                    dx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
                                    dy2*(phi[i][j-1][k] + phi[i][j+1][k]) +
                                    dz2*(phi[i][j][k-1] + phi[i][j][k+1])) / 
                                    (2.0 * (dx2+dy2+dz2));

                        // Apply Successive Over-Relaxation
                        phi[i][j][k] = phi[i][j][k] + 1.4 * (phi_new - phi[i][j][k]);
                    }
                }
            }

            // Check for convergence
            if iteration % 25 == 0 {
                let mut sum = 0.0;
                for i in 1..self.domain.nodes[0]-1 {
                    for j in 1..self.domain.nodes[1]-1 {
                        for k in 1..self.domain.nodes[2]-1 {
                            let r = -phi[i][j][k] * (2.0 * (dx2+dy2+dz2)) + 
                                rho[i][j][k] / PERMITTIVITY +
                                dx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
                                dy2*(phi[i][j-1][k] + phi[i][j+1][k]) +
                                dz2*(phi[i][j][k-1] + phi[i][j][k+1]);

                            sum += r*r;
                        }

                        let norm = (sum / (self.domain.nodes[0] * self.domain.nodes[1] * self.domain.nodes[2]) as f32).sqrt();
                        if norm < self.tolerance {
                            return true;
                        }
                    }
                }
            }

        }

        return false;
    }

    /// Computes the electric field $$-\nabla{\phi}$$
    pub fn compute_electric_field(&self) {
    }
}

