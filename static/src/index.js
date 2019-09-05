import Vue from 'vue';
import ElementUI from 'element-ui';
import 'element-ui/lib/theme-chalk/index.css';
import App from './pages/app.vue';



Vue.use(ElementUI);
new Vue({
    el:'#root',
    render:h=>h(App)
})
// //申明一个空div
// const root=document.createElement('div');
// document.body.appendChild(root);
//
// //把app.vue的内容挂载到空div上
// new Vue({
//     render:(h) =>h(App)
// }).$mount(root);