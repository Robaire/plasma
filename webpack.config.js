const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const WasmPackPlugin = require("@wasm-tool/wasm-pack-plugin");

const dist = path.resolve(__dirname, "dist");

module.exports = {
    mode: "development",
    entry: {
        index: "./js/index.js"
    },
    devServer: {
        contentBase: dist,
        writeToDisk: true
    },
    plugins: [
        new HtmlWebpackPlugin({
            title: "Plasma",
            template: "./static/index.html"
        }),
        new WasmPackPlugin({
            crateDirectory: path.resolve(__dirname, "."),
            outName: "pic"
        })
    ],
    output: {
        path: dist,
        filename: "[name].js",
        clean: true
    },
    experiments: {
        asyncWebAssembly: true
    }

};
