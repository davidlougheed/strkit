import HtmlWebpackPlugin from "html-webpack-plugin";
import path from "path";

export default {
    mode: "development",
    entry: "./src/index.tsx",
    module: {
        rules: [
            {
                test: /\.tsx?$/,
                use: "ts-loader",
                exclude: /node_modules/,
            },
            {
                test: /\.css$/i,
                use: ["style-loader", "css-loader"],
            },
        ],
    },
    resolve: {
        extensions: [".tsx", ".ts", ".js"],
    },
    output: {
        filename: "bundle.js",
        path: path.resolve(import.meta.dirname, "dist"),
    },
    plugins: [
        new HtmlWebpackPlugin({
            template: path.resolve(import.meta.dirname, "./src/template.html"),
            hash: true,
            publicPath: "/",
        }),
    ],
    devServer: {
        static: {
            directory: path.join(import.meta.dirname, "public"),
        },
        historyApiFallback: true,
    },
};
