const Vue = require("vue/dist/vue.js");
const { remote } = require("electron");
const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const process = require("process");

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
            in: "",
            output: "",
            mapping: "",
            term: "",
            type: "",
            species: "",
            interaction_confidence_score: 0.0,
            cluego_grouping: "",
            debug: false,
            cluego_p_value: 1.0,
        },
        stdout: "",
        stderr: "",
        running: false,
    },
    methods: {
        run: function() {
            let that = this;
            if(this.running) {
                return;
            }

            this.stdout = "";
            this.stderr = "";

            this.saveSession();

            const pine = spawn(path.join(__dirname, "/../../pine_2/pine_2.exe"), ["-i", this.input.in, "-o", path.join(this.input.output, "Merged_Input.csv"), "-t", this.input.type, "-s", this.input.species, "-m", this.input.mapping, "-v", path.join(this.input.output, "PINE.cys")]);

            pine.stdout.on("data", function(d) {
                that.stdout += d + "\n";
            });

            pine.stderr.on("data", function(d) {
                that.stderr += d + "\n";
            });

            pine.on("close", function(code) { 
                if(code == 0) {
                    that.stdout += "process completed successfully\n";
                } else {
                    that.stdout += "process failed\n";
                }
                that.running = false;
            });
        },
        inputFileChange: function(e, name) {
            if(!e.target.files[0].path) {
                return;
            }
            this.input[name] = e.target.files[0].path;
        },
        saveSession: function() {
            let session_file = this.getSessionFile();
            fs.writeFileSync(session_file, JSON.stringify(this.input));
        },
        loadSession: function() {
            let session_file = this.getSessionFile();
            if(!fs.existsSync(session_file)) {
                return;
            }
            let raw = fs.readFileSync(session_file, "utf8");
            let input = JSON.parse(raw);
            Vue.set(this, "input", input);
        },
        getSessionFile: function() {
            let app_dir = process.env.APPDATA || (process.platform == 'darwin' ? process.env.HOME + 'Library/Preferences' : process.env.HOME + "/.local/share");
            let full_dir = path.join(app_dir, "cedars-sinai-jve", "pine");
            if(!fs.existsSync(full_dir)) {
                fs.mkdirSync(full_dir, {recursive: true});
            }
            return path.join(full_dir, "session.json");
        },
    },
    mounted: function() {
        this.loadSession();
    },
});