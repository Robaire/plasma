use wasm_bindgen::prelude::*;
use std::ops::Index;

#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
extern {
    fn alert(s: &str);
}

#[wasm_bindgen]
pub fn init() {

    #[cfg(feature = "console_error_panic_hook")]
    std::panic::set_hook(Box::new(console_error_panic_hook::hook));

    #[cfg(feature = "wasm-logger")]
    wasm_logger::init(wasm_logger::Config::default());
}

#[wasm_bindgen]
pub fn greet(name: &str) {
    log::info!("Some info!");
    alert(&format!("Hello, {}!", name));
}

struct World {
    origin: [f32; 3], // Mesh Origin
    spacing: f32, // Cell Spacing
    nodes: u32, // Number of Nodes in each Axis
}

impl World {
    pub fn new(nodes: u32, spacing: f32) -> World {
        World {
            origin: [0.0, 0.0, 0.0],
            spacing: spacing,
            nodes: nodes
        }
    }
}

struct Field {
    data: Vec<Vec<Vec<f32>>>
}

impl Field {
    pub fn new(size: usize) -> Field {
        Field {
            data: vec![vec![vec![0.0; size]; size]; size]
        }
    }
}

impl Index<usize> for Field {
    type Output = Vec<Vec<f32>>;

    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}

