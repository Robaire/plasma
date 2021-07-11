use wasm_bindgen::prelude::*;

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
