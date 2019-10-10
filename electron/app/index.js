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
const CLUEGO_CONFIGURATION_BASE_NAME = "ClueGOConfiguration";

let vm = new Vue({
    el: "#app",
    data: {
        allowed_types: {"noFC": "No fold change", "singleFC": "Single fold change", "multiFC": "Multi fold change", "category": "Category"},
        species_map: {"Homo Sapiens": "human", "Mus Musculus": "mouse", "Rattus Norvegicus": "rat"},
        allowed_runs: ["string", "genemania", "both"],
        allowed_leading_terms: ["highest significance", "no. of genes per term", "percent of genes per term", "percent genes per term vs cluster"],
        allowed_visualize: ["biological process","subcellular location","molecular function","pathways","all"],
        allowed_grouping: ["global", "medium", "detailed"],

        input: {
            cytoscape_path: "",
            cluego_base_path: "",
            in: "",
            output: "",
            type: "",
            species: "",
            limit: 0,
            score: 0.4,
            significant: false,
            run: "both",
            fccutoff: 0.0,
            pvalcutoff: 1.0,
            leading_term: "no. of genes per term",
            visualize: "pathways",
            cluego_pval: 0.05,
            reference_path: "",
            grouping: "medium",
            debug: false,
        },
        automatic_input: { // for messaging the user that the paths were found automatically
            cytoscape_path: false,
            cluego_base_path: false,
        },
        stdout: "",
        stderr: "",
        error: "",
        running: false,
        tabs: TABS,
        current_tab: TABS.PATH_LOCATION,
        cluego_versions: null,
        cluego_picked_version: null,
    },
    methods: {
        run: function() {
            let that = this;
            if(!this.runnable()) {
                return;
            }
            if(this.running) {
                return;
            }
            this.running = true;
            this.current_tab = TABS.PROGRESS;

            this.stdout = "";
            this.stderr = "";

            this.saveSession();

            let args = [
                "--in", this.input.in,
                "--output", path.join(this.input.output, "Merged_Input.csv"),
                "--type", this.input.type,
                "--species", this.input.species,
                "--mapping", this.get_cluego_mapping(),
                "--save-session", path.join(this.input.output, "PINE.cys"),
                "--limit", this.limit,
                "--score", this.score,
                "--run", this.run,
                "--fccutoff", this.fccutoff,
                "--pvalcutoff", this.pvalcutoff,
                "--output_cluego", path.join(this.input.output, "output_cluego.txt"),
                "--leading-term", this.leading_term,
                "--visualize", this.visualize,
                "--cluego-pval", this.cluego_pval,
                "--reference-path", this.reference_path,
                "--grouping", this.grouping,
            ];
            if(this.significant) {
                args.push("--significant");
            }
            if(this.debug) {
                args.push("--debug");
            }
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
        runnable: function() {
            if(this.input.in && this.get_cluego_mapping() && this.input.output && this.input.type) {
                return true;
            }
            return false;
        },
        get_cluego_mapping: function() {
            if(!this.cluego_picked_version) {
                return null;
            }

            for(const mapping of this.cluego_picked_version.mapping_files) {
                if(this.species_map[mapping.species] === this.input.species) {
                    return mapping.fullpath;
                }
            }
            return null;
        },
        inputFileChange: function(e, name) {
            if(!e.target.files[0].path) {
                return;
            }

            let path = e.target.files[0].path;
            if(name === "cytoscape_path") {
                this.setCluegoBasePath(path);
                this.automatic_input.cytoscape_path = false;
            } else if(name === "cluego_base_path") {
                this.setCluegoBasePath(path);
                this.automatic_input.cluego_base_path = false;
            } else {
                this.input[name] = path;
            }
            this.refreshTab();
        },
        saveSession: function() {
            let session_file = this.getSessionFile();
            fs.writeFileSync(session_file, JSON.stringify(this.input));
        },
        loadSession: function() {
            const session_file = this.getSessionFile();
            if(!fs.existsSync(session_file)) {
                return;
            }
            const raw = fs.readFileSync(session_file, "utf8");
            const saved_input = JSON.parse(raw);
            for(const key in saved_input) {
                if(key in this.input) {
                    if(key === "cluego_base_path") {
                        this.setCluegoBasePath(saved_input[key]);
                    } else if(key === "cytoscape_path") {
                        this.setCytoscapePath(saved_input[key]);
                    } else {
                        Vue.set(this.input, key, saved_input[key]);
                    }
                }
            }
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
                        const cyto_path = path.join(fd, "Cytoscape.exe");
                        if(fs.existsSync(cyto_path)) {
                            this.setCytoscapePath(cyto_path);
                            this.automatic_input.cytoscape_path = true;
                            break;
                        }
                    }
                }
            }

            if(!this.input.cluego_base_path) {
                const cluego_base_path = path.join(os.homedir(), "ClueGOConfiguration");
                if(fs.existsSync(cluego_base_path)) {
                    this.setCluegoBasePath(cluego_base_path);
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
        setCluegoBasePath: function(cluego_path) {
            this.error = "";

            if(!fs.existsSync(cluego_path)) {
                this.error = "ClueGO configuration file path does not exist";
                return;
            }

            const path_split = cluego_path.split(path.sep);
            const path_index = path_split.indexOf(CLUEGO_CONFIGURATION_BASE_NAME);
            if(path_index === -1) {
                this.error = "ClueGO configuration directory must be named " + CLUEGO_CONFIGURATION_BASE_NAME;
                return;
            }

            cluego_path = path_split.slice(0, path_index + 1).join(path.sep);
            if(!fs.statSync(cluego_path).isDirectory()) {
                this.error = CLUEGO_CONFIGURATION_BASE_NAME + " must be a directory";
                return;
            }

            // find all valid versions of cluego to use
            let version_dirs = [];
            let sub_dirs = fs.readdirSync(cluego_path);
            for(const version of sub_dirs) {
                const sd = path.join(cluego_path, version);
                if(!fs.statSync(sd).isDirectory()) {
                    continue;
                }
                let source_file_dir = path.join(sd, "ClueGOSourceFiles");
                if(!fs.existsSync(source_file_dir) || !fs.statSync(source_file_dir).isDirectory()) {
                    /* source files must exist and be a directory */
                    continue;
                }

                let mapping_files = []
                let source_files_list = fs.readdirSync(source_file_dir);
                for(const sf of source_files_list) {
                    if(!sf.startsWith("Organism_")) {
                        continue;
                    }
                    let species_name = sf.replace("Organism_", "");
                    if(!(species_name in this.species_map)) {
                        continue;
                    }
                    const sf_full = path.join(source_file_dir, sf);
                    if(!fs.statSync(sf_full).isDirectory()) {
                        continue;
                    }

                    const accession_files = fs.readdirSync(sf_full);
                    let accession_file = null;
                    for(const af of accession_files) {
                        if(af.includes("gene2accession") && af.endsWith(".txt.gz") && !fs.statSync(path.join(sf_full, af)).isDirectory()) {
                            accession_file = af;
                            break;
                        }
                    }
                    if(accession_file !== null) {
                        mapping_files.push({
                            species: species_name,
                            filename: accession_file,
                            fullpath: path.join(sf_full, accession_file),
                        });
                    }
                }

                if(mapping_files.length > 0) {
                    version_dirs.push({
                        version: version,
                        mapping_files: mapping_files,
                        latest: false,
                    });
                }
            }

            if(version_dirs.length > 0) {
                this.input.cluego_base_path = cluego_path;
                this.cluego_versions = version_dirs;
                this.cluego_versions.sort(function(a, b) {
                    return -1 * a.localeCompare(b);
                });
                this.cluego_picked_version = this.cluego_versions[0];
                this.cluego_picked_version.latest = true;
            }
        },
        setCytoscapePath: function(path) {
            if(fs.existsSync(path) && fs.statSync(path).isFile()) {
                this.input.cytoscape_path = path;
            }
        },
    },
    mounted: function() {
        this.loadSession();
        this.searchForPaths();
    },
    filters: {
        filename: function(v) {
            return path.basename(v);
        }
    }
});