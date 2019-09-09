<template>
    <el-container class="is-vertical">
        <el-row>
            <el-col :span="8">
                <el-card class="box-card">
                    <div slot="header" class="clearfix">
                        <span>基础信息</span>
                    </div>
                    <div>
                        <el-form ref="baseInfoForm" label-position="left" label-width="80px" :model="baseInfoForm">
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
                                <el-select v-model="baseInfoForm.WellType" placeholder="井类型">
                                    <el-option label="直井" value="vertical"></el-option>
                                </el-select>
                            </el-form-item>
                            <el-form-item label="油气类型">
                                <el-select v-model="baseInfoForm.OilGasType" placeholder="油气类型">
                                    <el-option label="油" value="oil"></el-option>
                                </el-select>
                            </el-form-item>

                        </el-form>

                    </div>
                </el-card>
            </el-col>

            <el-col :span="8">
                <el-card class="box-card">
                    <div slot="header" class="clearfix">
                        <span>基础信息</span>
                    </div>
                    <div>
                        <el-form ref="baseInfoForm" label-position="left" label-width="80px" :model="baseInfoForm">
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
                                <el-select v-model="baseInfoForm.WellType" placeholder="井类型">
                                    <el-option label="直井" value="vertical"></el-option>
                                </el-select>
                            </el-form-item>
                            <el-form-item label="油气类型">
                                <el-select v-model="baseInfoForm.OilGasType" placeholder="油气类型">
                                    <el-option label="油" value="oil"></el-option>
                                </el-select>
                            </el-form-item>

                        </el-form>

                    </div>
                </el-card>
            </el-col>

            <el-col :span="8">
                <el-card class="box-card">
                    <div slot="header" class="clearfix">
                        <span>基础信息</span>
                    </div>
                    <div>
                        <el-form ref="baseInfoForm" label-position="left" label-width="80px" :model="baseInfoForm">
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
                                <el-select v-model="baseInfoForm.WellType" placeholder="井类型">
                                    <el-option label="直井" value="vertical"></el-option>
                                </el-select>
                            </el-form-item>
                            <el-form-item label="油气类型">
                                <el-select v-model="baseInfoForm.OilGasType" placeholder="油气类型">
                                    <el-option label="油" value="oil"></el-option>
                                </el-select>
                            </el-form-item>

                        </el-form>

                    </div>
                </el-card>
            </el-col>
        </el-row>
        <el-row>
            <div style="margin:20px;">
            </div>
        </el-row>
        <el-row>
            <div style="margin:0 auto;width:200px;">
                <el-button type="primary" @click="onSubmit">立即创建</el-button>
                <el-button @click="onReset">重置</el-button>
            </div>
        </el-row>
    </el-container>
</template>
<script>
    export default {
        data() {
            return {

                baseInfoForm: {
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
                ,
                onReset() {
                    console.log('reset!');
                }
            }
        ,
        computed: {}

    }

</script>
<style scoped>

</style>