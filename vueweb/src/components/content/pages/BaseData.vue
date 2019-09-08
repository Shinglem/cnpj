<template>

    <el-form ref="baseInfoForm" :model="baseInfoForm" labl-width="80px">
        <el-form-item label="油田名称">
            <el-input v-model="baseInfoForm.FieldName"></el-input>
        </el-form-item>
        <el-form-item label="井名">
            <el-input v-model="baseInfoForm.WellName"></el-input>
        </el-form-item>
        <el-form-item label="产层">
            <el-input v-model="baseInfoForm.Zone"></el-input>
        </el-form-item>
        <el-form-item label="井类型">
            <el-select v-model="form.WellType" placeholder="井类型">
                <el-option label="直井" value="vertical"></el-option>
            </el-select>
        </el-form-item>
        <el-form-item label="油气类型">
            <el-select v-model="form.OilGasType" placeholder="油气类型">
                <el-option label="油" value="oil"></el-option>
            </el-select>
        </el-form-item>
        <el-form-item>
            <el-button type="primary" @click="onSubmit">立即创建</el-button>
        </el-form-item>

    </el-form>


</template>
<script>
    export default {
        data() {
            return {

                form: {
                    FieldName: '',
                    WellName: '',
                    Zone: '',
                    WellType: '',
                    OilGasType: ''

                }


            }
        }
        ,

        methods:
            {
                onSubmit() {
                    console.log('submit!');
                }
            }
        ,
        computed: {
            activeName() {
                return this.$store.state.test
            }
        }

    }

</script>
<style scoped>

</style>