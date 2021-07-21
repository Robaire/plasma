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

    let domain = Domain::new([100, 100, 100], [0.1, 0.1, 0.1]);

    let mut field = Field::new_cube(100);

    log::info!("Element: {}", field[0][0][0]);

    let b = &mut field[0][0][0];
    *b = 3.14;

    log::info!("Element: {}", field[0][0][0]);
}

/// Represents the computational domain of the simulation
struct Domain {
    nodes: [u32; 3], // The number of nodes in each axis
    spacing: [f32; 3], // Space between nodes in x, y, and z
}

impl Domain {

    /// Create a new domain
    pub fn new(nodes: [u32; 3], size: [f32; 3]) -> Domain {
        Domain {
            nodes: nodes,
            spacing: [
                size[0] / nodes[0] as f32,
                size[1] / nodes[1] as f32,
                size[2] / nodes[2] as f32
            ]
        }
    }
}

/// 3D field of floats
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

