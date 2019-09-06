<template>
    <el-container style="height: 500px; border: 1px solid #eee">

        <el-header style="text-align: center; font-size: 20px">
            <Head/>
        </el-header>

        <el-container>
            <el-aside width="200px" style="background-color: rgb(238, 241, 246)">
                <Left/>
            </el-aside>


            <el-main>
                <Main/>
            </el-main>
        </el-container>
    </el-container>


</template>

<style>
    .el-header {
        background-color: #B3C0D1;
        color: #333;
        line-height: 60px;
    }

    .el-aside {
        color: #333;
    }
</style>

<script>

    import Head from './components/frame/Head.vue'
    import Left from './components/frame/Left.vue'
    import Main from './components/content/Main.vue'



    export default {

        components: {
            Head, Left, Main
        }
    };
</script>