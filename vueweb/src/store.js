import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
  state: {
    test : 'first'
  },
  mutations: {
    setTest(state , test){
      state.test = test
    }
  },
  actions: {

  }
})
