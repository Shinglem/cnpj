const path = require('path');
const {VueLoaderPlugin} = require("vue-loader");
const webpack = require('webpack');
module.exports = {
    entry: {
        app: './src/index.js'
    }, // 项目的入口文件，webpack会从main.js开始，把所有依赖的js都加载打包
    output: {
        filename: '[name].build.js', // 打包后的文件名
        path: path.resolve(__dirname, './dist'), // 项目的打包文件路径
        //publicPath: '/dist/' // 通过devServer访问路径
    },
    module: {
        rules: [
            {
                test: /\.vue$/,
                use: ['vue-loader'] // +++
            },
            {
                test: /\.css$/,
                use: ['vue-style-loader', 'css-loader'] // +++
            },

            {
                test: /\.(eot|woff|ttf)$/,
                loader: 'file-loader'
            },

            {
                test: /\.html$/,
                loader: 'html-loader'
            }
        ]
    },
    plugins: [
        new VueLoaderPlugin() // +++
    ]
};