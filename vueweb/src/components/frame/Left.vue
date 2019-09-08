<template>
    <div>
        <el-menu>
            <el-submenu index="1">
                <template slot="title"><i class="el-icon-menu"></i>地质评价</template>
                <el-menu-item-group>
                    <el-menu-item index="1-1">测井数据录入</el-menu-item>
                    <el-menu-item index="1-2">测井曲线剖面</el-menu-item>

                </el-menu-item-group>

            </el-submenu>
            <el-submenu index="2">
                <template slot="title"><i class="el-icon-menu"></i>非稳态产能-EUR评价</template>
                <el-menu-item-group>
                    <el-menu-item index="2-1">压力数据录入</el-menu-item>
                    <el-menu-item index="2-2">产量数据录入</el-menu-item>
                    <el-menu-item index="2-3">数据预处理</el-menu-item>
                    <el-menu-item index="2-4">井口压力折算</el-menu-item>
                    <el-menu-item index="2-5">解释段选择</el-menu-item>
                    <el-menu-item index="2-6">参数解释拟合</el-menu-item>
                    <el-menu-item index="2-7">自动拟合</el-menu-item>
                    <el-menu-item index="2-8">生产预测</el-menu-item>
                </el-menu-item-group>

            </el-submenu>
            <el-submenu index="3">
                <template slot="title"><i class="el-icon-menu"></i>稳态产能评价</template>
                <el-menu-item-group>
                    <el-menu-item index="3-1">地质-工程参数录入</el-menu-item>
                    <el-menu-item index="3-2">产能计算</el-menu-item>
                    <el-menu-item index="3-3">生产预测
                    </el-menu-item>
                </el-menu-item-group>
            </el-submenu>

            <el-submenu index="4">
                <template slot="title"><i class="el-icon-menu"></i>产能试井评价</template>
                <el-menu-item-group>
                    <el-menu-item index="4-1">导入产能试井数据</el-menu-item>
                    <el-menu-item index="4-2">解释校正</el-menu-item>
                    <el-menu-item index="4-3">生产预测</el-menu-item>
                </el-menu-item-group>
            </el-submenu>

            <el-submenu index="5">
                <template slot="title"><i class="el-icon-menu"></i>EUR评价</template>
                <el-menu-item-group>
                    <el-menu-item index="5-1">物质平衡方法</el-menu-item>
                    <el-menu-item index="5-2">流动物质平衡</el-menu-item>
                    <el-menu-item index="5-3">递减分析</el-menu-item>
                    <el-menu-item index="5-3">经验公式</el-menu-item>
                </el-menu-item-group>
            </el-submenu>
        </el-menu>
    </div>
</template>

<script>
    export default {
        name: "Left"
    }
</script>

<style scoped>

</style>