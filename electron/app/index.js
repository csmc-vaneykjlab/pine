let Vue = require("vue/dist/vue.js");
const remote = require("electron").remote;

document.addEventListener("keydown", function(e) {
    if(e.keyCode === 123) {
        let window = remote.getCurrentWindow();
        window.toggleDevTools();
    }
});

let vm = new Vue({
    el: "#app",
    data: {
        allowed_terms: ["pathways", "biological process", "subcellular location", "molecular function", "all"],
        allowed_types: ["noFC", "singleFC", "multiFC", "category"],
        allowed_species: ["human", "mouse", "rat"],
        allowed_cluego_grouping: ["global", "medium", "detailed"],
        input: {
            input: "",
            output: "",
            cluego_mapping: "",
            term: "",
            type: "",
            species: "",
            interaction_confidence_score: 0.0,
            cluego_grouping: "",
            debug: false,
            cluego_p_value: 1.0,
        }
    },
    methods: {

    }
});