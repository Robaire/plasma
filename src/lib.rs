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

/// Computes the Debye length in UNIT?
///
/// * `temp` - Species temperature in eV
/// * `density` - Species density in # / m^3
fn debye_length(temp: f32, density: f32) -> f32 {

    // 1 eV = (q/k_b) 
    // 1 eV ~= 11604.519 K
    // Since temp is expresssed in eV, the boltzmann constant cancels out

    let permittivity: f32 = 55.263_494_06 * 10e6; // Multiple by 10e6 to convert \frac{e^2}{GeV fm} -> \frac{e^2}{eV m}
    // let boltzmann: f32 = 8.617_333_262_145e-5; // \frac{eV}{K}
    // let charge: f32 = -1.0; // {e}

    // Due to some unit conversions we can simplify the calculation to this
    return ((permittivity * temp)/(density)).sqrt();
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

    pub fn new(domain: &mut Domain, iterations: usize, tolerance: f32) -> PotentialSolver {
        PotentialSolver {
            domain,
            iterations,
            tolerance
        }
    }

    pub fn solve(&mut self) -> bool {

        let ep: f32 = 1.0; // Permittivity of free space

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
                        let phi_new = ((rho[i][j][k]/ep) + 
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
                                rho[i][j][k] / ep +
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

