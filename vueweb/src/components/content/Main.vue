<template>
    <div>

        <el-tabs activeName="first" @tab-click="handleClick">
            <el-tab-pane label="基础数据" name="first">
                <BaseData/>
            </el-tab-pane>
            <!--            <el-tab-pane label="配置管理" name="second">配置管理</el-tab-pane>-->
            <!--            <el-tab-pane label="角色管理" name="third">角色管理</el-tab-pane>-->
            <!--            <el-tab-pane label="定时任务补偿" name="fourth">定时任务补偿</el-tab-pane>-->
        </el-tabs>


    </div>
</template>
<script>
    import BaseData from './pages/BaseData.vue'

    export default {

        methods: {
            handleClick(tab, event) {
                console.log(tab, event);
            }
        },
        computed: {


            // activeName: {
            //     // getter
            //     get: function() {
            //         return this.$store.state.test
            //     }
            // ,
            //     // setter
            //     set: function(newValue) {
            //         this.$store.state.test = newValue
            //     }
            //
            // }
        }
        ,
        components: {
            BaseData
        }
    };
</script>
<style scoped>

</style>