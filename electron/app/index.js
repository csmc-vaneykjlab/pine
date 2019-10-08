const Vue = require("vue/dist/vue.js");
const { remote } = require("electron");
const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const process = require("process");
const os = require("os");

document.addEventListener("keydown", function(e) {
    if(e.keyCode === 123) {
        let window = remote.getCurrentWindow();
        window.toggleDevTools();
    }
});

const TABS = {
    PATH_LOCATION: "path-location",
    INPUT: "input",
    PROGRESS: "progress",
    PATHWAY_SELECTION: "pathway-selection",
    COMPLETE: "complete",
};

let vm = new Vue({
    el: "#app",
    data: {
        allowed_terms: ["pathways", "biological process", "subcellular location", "molecular function", "all"],
        allowed_types: ["noFC", "singleFC", "multiFC", "category"],
        allowed_species: ["human", "mouse", "rat"],
        allowed_cluego_grouping: ["global", "medium", "detailed"],
        input: {
            cytoscape_path: "",
            cluego_base_path: "",
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
        automatic_input: { // for messaging the user that the paths were found automatically
            cytoscape_path: false,
            cluego_base_path: false,
        },
        stdout: "",
        stderr: "",
        running: false,
        tabs: TABS,
        current_tab: TABS.PATH_LOCATION,
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

            let args = ["-i", this.input.in, "-o", path.join(this.input.output, "Merged_Input.csv"), "-t", this.input.type, "-s", this.input.species, "-m", this.input.mapping, "-v", path.join(this.input.output, "PINE.cys")];
            if(process.env.NODE_ENV === "dev") {
                var pine = spawn("C:/Users/GoJ1/AppData/Local/Programs/Python/Python37/python.exe", [path.join(__dirname, "/../../pine_2.py")].concat(args));
            } else {
                var pine = spawn(path.join(__dirname, "/../../pine_2/pine_2.exe"), args);
            }

            pine.stdout.on("data", function(d) {
                that.stdout += d + "\n";
            });

            pine.stderr.on("data", function(d) {
                that.stderr += d + "\n";
            });

            pine.on("close", function(code) { 
                if(code === 0) {
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
            if(name === "cytoscape_path") {
                this.automatic_input.cytoscape_path = false;
            } else if (name === "cluego_base_path") {
                this.automatic_input.cluego_base_path = false;
            }
            this.refreshTab();
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
        searchForPaths: function() {
            /* only support for windows */
            if(os.platform() !== "win32") {
                return;
            }

            if(!this.input.cytoscape_path) {
                const programs_path = "C:/Program Files";
                let found_dirs = [];
                if(fs.existsSync(programs_path)) {
                    const files = fs.readdirSync(programs_path);
                    for(const file of files) {
                        if(!file.startsWith("Cytoscape_v3")) {
                            continue;
                        }
                        let full_path = path.join(programs_path, file);
                        if(!fs.statSync(full_path).isDirectory()) {
                            continue;
                        }
                        found_dirs.push(path.join(full_path));
                    }
                }
                if(found_dirs.length > 0) {
                    /* reverse sort */
                    found_dirs.sort(function(a, b) {
                        return -1 * a.localeCompare(b);
                    });

                    for(const fd of found_dirs) {
                        const cyto_path = path.join(fd, "cytoscape.exe");
                        if(fs.existsSync(cyto_path)) {
                            this.input.cytoscape_path = cyto_path;
                            this.automatic_input.cytoscape_path = true;
                            break;
                        }
                    }
                }
            }

            if(!this.input.cluego_base_path) {
                const cluego_base_path = path.join(os.homedir(), "ClueGOConfiguration");
                if(fs.existsSync(cluego_base_path)) {
                    //this.setCluegoBasePath(cluego_base_path);
                    this.input.cluego_base_path = cluego_base_path;
                    this.automatic_input.cluego_base_path = true;
                }
            }

            this.refreshTab();
        },
        refreshTab: function() {
            if(this.current_tab === TABS.PATH_LOCATION && this.input.cytoscape_path && this.input.cluego_base_path) {
                this.current_tab = TABS.INPUT;
            }
        },
        setCluegoBasePath: function(path) {
        },
    },
    mounted: function() {
        this.loadSession();
        this.searchForPaths();
    },
});