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
    import Head from './side/head.vue'
    import Left from './side/left.vue'
    import Main from './side/main.vue'

    export default {
        // data() {
        //     const item = {
        //         date: '2016-05-02',
        //         name: '王小虎',
        //         address: '上海市普陀区金沙江路 1518 弄'
        //     };
        //     return {
        //         tableData: Array(20).fill(item)
        //     }
        // },
        components: {
            Head , Left , Main
        }
    }
</script>

<!--<script>-->
<!--    export default {-->
<!--        data() {-->
<!--            const item = {-->
<!--                date: '2016-05-02',-->
<!--                name: '王小虎',-->
<!--                address: '上海市普陀区金沙江路 1518 弄'-->
<!--            };-->
<!--            return {-->
<!--                tableData: Array(20).fill(item)-->
<!--            }-->
<!--        }-->
<!--    };-->
<!--</script>-->