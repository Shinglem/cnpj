<template>
    <div>
        <el-menu class="el-menu-demo" mode="horizontal">
            <el-submenu index="1">
                <template slot="title"><i class="el-icon-menu"></i>项目</template>
                <el-menu-item index="1-1">新建</el-menu-item>
                <el-menu-item index="1-2">打开</el-menu-item>
                <el-menu-item index="1-3">保存</el-menu-item>
                <el-menu-item index="1-4">属性</el-menu-item>
                <el-menu-item index="1-5">退出</el-menu-item>
            </el-submenu>

            <el-submenu index="2">
                <template slot="title"><i class="el-icon-menu"></i>视图</template>
                <el-menu-item index="2-1">选择</el-menu-item>
                <el-menu-item index="2-2">放大</el-menu-item>
                <el-menu-item index="2-3">缩小</el-menu-item>
                <el-menu-item index="2-4">复位</el-menu-item>
            </el-submenu>

        </el-menu>
        <div class="line"></div>

    </div>
</template>

<script>
    export default {
        name: "Options"
    }
</script>

<style scoped>

</style>