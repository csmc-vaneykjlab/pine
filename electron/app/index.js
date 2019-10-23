const Vue = require("vue/dist/vue.js");
const { remote, net } = require("electron");
const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const process = require("process");
const os = require("os");
const readline = require("readline");
const shell = require("electron").shell;

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
};
const CLUEGO_CONFIGURATION_BASE_NAME = "ClueGOConfiguration";

let vm = new Vue({
    el: "#app",
    data: {
        allowed_types: {
            "noFC": {"text": "No fold change", "er": false}, // er - extra required - extra fields need to be shown and submitted by the user
            "singleFC": {"text": "Single fold change", "er": false},
            "multiFC": {"text": "Multi fold change", "er": false},
            "category": {"text": "Category", "er": false},
            "singleFC-site": {"text": "Single fold change site", "er": true},
            "multiFC-site": {"text": "Multi fold change site", "er": true},
        },
        species_map: {"homo sapiens": "human", "mus musculus": "mouse", "rattus norvegicus": "rat"},
        allowed_runs: ["string", "genemania", "both"],
        allowed_visualize: ["biological process","subcellular location","molecular function","pathways","all"],
        allowed_grouping: ["global", "medium", "detailed"],
        allowed_enzymes: {
            "trypsin": "Trypsin",
            "trypsin_p": "Trypsin P",
            "lys_n": "Lysine N-terminus",
            "lys_c": "Lysine C-terminus",
            "asp_n": "Asparagine N-terminus",
            "arg_n": "Arginine N-terminus",
            "chymotrypsin": "Chymotrypsin",
        },

        input: {
            cytoscape_path: "",
            cluego_base_path: "",
            in: "",
            output: "",
            type: null,
            species: null,
            limit: null,
            score: null,
            significant: null,
            run: null,
            fccutoff: null,
            pvalcutoff: null,
            visualize: null,
            cluego_pval: null,
            reference_path: null,
            grouping: null,
            enzyme: null,
            fasta_file: null,
            mods: null,
        },
        automatic_input: { // for messaging the user that the paths were found automatically
            cytoscape_path: false,
            cluego_base_path: false,
        },
        session_files: {
            cytoscape_session: null,
            cluego_pathways: null,
        },
        stdout: "",
        stderr: "",
        error: "",
        running: false,
        tabs: TABS,
        current_tab: TABS.PATH_LOCATION,
        cluego_versions: null,
        cluego_picked_version: null,
        cluego_pathways: {
            data: null,
            header: null,
            query: null,
            page: null,
            per_page: null,
            sources: null,
            sort: null,
        },
        show_config: false,
        show_about: false,
    },
    methods: {
        run: async function(args) {
            let that = this;
            if(!this.runnable() || this.running) {
                return;
            }
            this.running = true;

            /* warn user about cytoscape running and give them a chance to stop */
            if(await this.is_cytoscape_running()) {
                const res = remote.dialog.showMessageBoxSync({
                    "type": "warning",
                    "buttons": [
                        "Continue anyways",
                        "Cancel (recommended)",
                    ],
                    "defaultId": 1,
                    "title": "Cytoscape open",
                    "message": "It is recommended that you close Cytoscape before continuing to avoid data loss",
                });
                if(res !== 0) {
                    this.running = false;
                    return;
                }
            }

            this.switchTab(TABS.PROGRESS);
            this.stdout = "Starting PINE analysis...\n\n";
            this.stderr = "";
            this.reset_cluego_pathways();

            this.saveSession();

            if(process.env.NODE_ENV === "dev") {
                let args1 = [path.join(__dirname, "/../../changes_to_pine_final.py")].concat(args);
                var pine = spawn("C:/Users/GoJ1/AppData/Local/Programs/Python/Python37/python.exe", args1);
            } else {
                var pine = spawn(path.join(__dirname, "/../../pine_2/pine_2.exe"), args);
            }

            pine.stdout.on("data", function(d) {
                if(typeof d !== "string") {
                    d = d.toString("utf8");
                }
                const d_split = d.split("\n");
                for(let ds of d_split) {
                    ds = ds.replace(/^\s+|\s+$/g, ""); // strip whitespace
                    if(ds.startsWith("COMMAND")) {
                        if(ds.startsWith("COMMAND FILE-PATHWAYS ")) {
                            that.session_files.cluego_pathways = ds.replace(/^COMMAND FILE-PATHWAYS /, "");
                        } else if(ds.startsWith("COMMAND FILE-SESSION ")) {
                            that.session_files.cytoscape_session = ds.replace(/^COMMAND FILE-SESSION /, "");
                        }
                    } else {
                        that.stdout += ds + "\n";
                    }
                }
            });

            pine.stderr.on("data", function(d) {
                that.stderr += d + "\n";
            });

            pine.on("close", function(code) { 
                that.running = false;
                if(code === 0) {
                    that.stdout += "run completed successfully\n";
                    that.read_cluego_pathways();
                    that.switchTab(TABS.PATHWAY_SELECTION);
                } else {
                    that.stdout += "run failed\n";
                    let res = remote.dialog.showMessageBoxSync({
                        "type": "warning",
                        "buttons": [
                            "Restart",
                            "Cancel",
                        ],
                        "defaultID": 1,
                        "title": "Run failed",
                        "message": that.stderr,
                        "noLink": true,
                    });
                    if(res === 0) { // restart
                        if(that.session_exists()) {
                            that.run_with_cluego_subset();
                        } else {
                            that.run_full();
                        }
                    } else { // cancel
                        if(that.session_exists()) {
                            that.switchTab(TABS.PATHWAY_SELECTION);
                        } else {
                            that.switchTab(TABS.INPUT);
                        }
                    }
                }
            });
        },
        run_full: function() {
            let args = this.pine_args();
            this.run(args);
        },
        run_with_cluego_subset: function() {
            if(!this.runnable_with_cluego_subset()) {
                return;
            }

            /* create the filtered cluego file */
            const filtered_file_name = this.session_files.cluego_pathways.replace(/\.txt$/, "") + "-filtered.txt";
            let file_data = [];
            file_data.push(this.cluego_pathways.header.join("\t"));
            for(const pathway of this.cluego_pathways.data) {
                if(pathway.selected) {
                    file_data.push(pathway.line);
                }
            }
            fs.writeFileSync(filtered_file_name, file_data.join("\n"), "utf-8");

            let args = this.pine_args();
            args.push("-a");
            args.push(filtered_file_name);
            args.push("-c");
            args.push(this.session_files.cytoscape_session);
            this.run(args);
        },
        reset_session_files: function() {
            this.session_files.cytoscape_session = null;
            this.session_files.cluego_pathways = null;
        },
        session_exists: function() {
            if(this.session_files.cytoscape_session || this.session_files.cluego_pathways) {
                return true;
            }
            return false;
        },
        is_cytoscape_running: async function() {
            /* check if listening on port 1234 */
            let running = await new Promise(function(resolve) {
                const request = remote.net.request("http://localhost:1234/v1/version");
                request.on("response", function() {
                    resolve(true);
                });
                request.on("error", function() {
                    resolve(false);
                });
                request.on("close", function() {
                    /* catch all, just return false */
                    resolve(false);
                });
                request.end();
            });
            return running;
        },
        pine_args: function() {
            let args = [
                "--in", this.input.in,
                "--output", this.input.output,
                "--type", this.input.type,
                "--species", this.input.species,
                "--mapping", this.get_cluego_mapping(),
                "--limit", this.input.limit,
                "--score", this.input.score,
                "--run", this.input.run,
                "--fccutoff", this.input.fccutoff,
                "--pvalcutoff", this.input.pvalcutoff,
                "--visualize", this.input.visualize,
                "--cluego-pval", this.input.cluego_pval,
                "--reference-path", this.input.reference_path,
                "--grouping", this.input.grouping,
                "--cytoscape-executable", this.input.cytoscape_path,
                "--enzyme", this.input.enzyme,
                "--fasta-file", this.input.fasta_file,
                "--mods", this.input.mods,
                "--gui",
            ];
            if(this.significant) {
                args.push("--significant");
            }

            return args;
        },
        runnable: function() {
            if(this.input.in && this.get_cluego_mapping() && this.input.output && this.input.type && !this.running) {
                if(this.is_extra_options_required()) {
                    if(this.input.mods && this.input.fasta_file && this.input.enzyme) {
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    return true;
                }
            }
            return false;
        },
        runnable_with_cluego_subset: function() {
            return this.runnable() && this.session_exists() && this.cluego_pathways.data.filter((x) => x.selected === true).length > 0;
        },
        reset_cluego_pathways: function() {
            this.cluego_pathways.data = [];
            this.cluego_pathways.header = [];
            this.cluego_pathways.query = "";
            this.cluego_pathways.page = 1;
            this.cluego_pathways.per_page = 5;
            this.cluego_pathways.sources = [];
            this.cluego_pathways.sort = null;
        },
        set_input_defaults: function() {
            this.input.type = "";
            this.input.species = "";
            this.input.limit = 0;
            this.input.score = 0.4;
            this.input.significant = false;
            this.input.run = "both";
            this.input.fccutoff = 0.0;
            this.input.pvalcutoff = 1.0;
            this.input.visualize = "pathways";
            this.input.cluego_pval = 0.05;
            this.input.reference_path = "";
            this.input.grouping = "medium";
            this.input.enzyme = "trypsin";
            this.input.fasta_file = "";
            this.input.mods = "S,T,Y";
        },
        read_cluego_pathways: function() {
            var that = this;

            this.reset_cluego_pathways();

            if(!this.session_files.cluego_pathways || !fs.existsSync(this.session_files.cluego_pathways) || !fs.statSync(this.session_files.cluego_pathways).isFile()) {
                return;
            }

            let counter = 0;
            let cluego_pathways = [];
            let line_reader =  readline.createInterface({
                input: fs.createReadStream(this.session_files.cluego_pathways),
            });
            let sources = new Set();
            let source_types = [
                {"phrase": "biologicalprocess", "name": "Biological Process"},
                {"phrase": "cellularcomponent", "name": "Cellular Component"},
                {"phrase": "molecularfunction", "name": "Molecular Function"},
                {"phrase": "human-diseases", "name": "Pathway (Human Diseases)"},
                {"phrase": "kegg", "name": "Pathway (KEGG)"},
                {"phrase": "wikipathways", "name": "Pathway (Wiki Pathways)"},
                {"phrase": "pathways", "name": "Pathway (Other)"},
            ];
            line_reader.on("line", function(line) {
                line = line.replace(/^\s+|\s+$/g, '');
                let fields = line.split("\t");

                if(counter === 0) {
                    that.cluego_pathways.header = fields;
                } else {
                    if(fields.length !== that.cluego_pathways.header.length) {
                        return; // invalid line
                    }
                    let record = {"data": {}, "selected": false, "line": line};
                    for(let i = 0; i < that.cluego_pathways.header.length; i++) {
                        record["data"][that.cluego_pathways.header[i]] = fields[i];
                    }
                    // don't split genes unless needed
                    // record["data"]["Associated Genes Found"] = record["data"]["Associated Genes Found"].replace(/^\[|\]$/g, '').split("\t");

                    /* get the ontology source */
                    let record_source = record["data"]["Ontology Source"].toLowerCase();
                    for(const st of source_types) {
                        //record_source.includes();
                    }
                }
                counter += 1;
            });

            this.cluego_pathways.data = cluego_pathways;
        },
        get_cluego_mapping: function() {
            if(!this.cluego_picked_version) {
                return null;
            }

            for(const mapping of this.cluego_picked_version.mapping_files) {
                if(this.species_map[mapping.species.toLowerCase()] === this.input.species) {
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
                this.setCytoscapePath(path);
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
                this.switchTab(TABS.INPUT);
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
                    if(!(species_name.toLowerCase() in this.species_map)) {
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
                    return -1 * a.version.localeCompare(b.version);
                });
                this.cluego_picked_version = this.cluego_versions[0];
                this.cluego_picked_version.latest = true;
            }
        },
        setCytoscapePath: function(cy_path) {
            if(!fs.existsSync(cy_path) || !fs.statSync(cy_path).isFile()) {
                return;
            }
            this.input.cytoscape_path = cy_path;
        },
        openGeneratedFile: function() {
            if(!fs.existsSync(this.session_files.cytoscape_session) || !fs.statSync(this.session_files.cytoscape_session).isFile()) {
                return;
            }
            shell.openItem(this.session_files.cytoscape_session);
        },
        tabReachable: function(tab_name) {
            switch(tab_name) {
                case TABS.PATH_LOCATION:
                    return false;
                case TABS.INPUT:
                    return true;
                case TABS.PROGRESS:
                    return this.running || this.session_exists();
                case TABS.PATHWAY_SELECTION:
                    return !this.running && this.session_exists();
            }
            return false;
        },
        switchTab: function(tab_name) {
            if(this.tabReachable(tab_name)) {
                this.current_tab = tab_name;
            }
        },
        run_disabled: function() {
            return this.running || this.session_files.cytoscape_session;
        },
        change_cluego_pathways_page: function(next_page) {
            if(next_page < 1) {
                return;
            }
            if(next_page > this.cluego_pathways_filtered.n_pages) {
                return;
            }
            this.cluego_pathways.page = next_page;
        },
        is_extra_options_required: function() {
            if(!(this.input.type in this.allowed_types)) {
                return false;
            }
            return this.allowed_types[this.input.type].er;
        },
    },
    mounted: function() {
        this.reset_cluego_pathways();
        this.set_input_defaults();
        this.loadSession();
        this.searchForPaths();

        this.session_files.cytoscape_session = "C:\\Users\\GoJ1\\Documents\\cytoscape1\\out\\20191022T145051-pine.cys"; // DEBUG:
        this.session_files.cluego_pathways = "C:\\Users\\GoJ1\\Documents\\cytoscape1\\out\\20191022T145051-cluego-pathways.txt"; // DEBUG:
        this.read_cluego_pathways();
        this.switchTab(TABS.PATHWAY_SELECTION); // DEBUG:
    },
    filters: {
        filename: function(v) {
            let basename = path.basename(v);
            if(basename.length > 10) {
                return basename.slice(0, 7) + "...";
            }
            return basename;
        },
        roundString: function(v, places) {
            const n = Number(v);
            if(isNaN(n)) {
                return v;
            } else {
                return n.toPrecision(places);
            }
        }
    },
    computed: {
        cluego_pathways_filtered: function() {
            var that = this;

            const filtered = this.cluego_pathways.data.filter(function(pathway) {
                if(that.cluego_pathways.query === "") {
                    return true;
                }

                let lower = that.cluego_pathways.query.toLowerCase();
                return pathway.data['GOTerm'].toLowerCase().includes(lower) || pathway.data['Ontology Source'].toLowerCase().includes(lower);
            });

            const full_length = filtered.length;
            const n_pages = Math.ceil(full_length / this.cluego_pathways.per_page);

            const start_idx = (this.cluego_pathways.page - 1) * this.cluego_pathways.per_page;
            const end_idx = start_idx + this.cluego_pathways.per_page;

            const data = filtered.slice(start_idx, end_idx);

            return {
                data: data,
                full_length: full_length,
                n_pages: n_pages,
            };
        },
        configuration_message: function() {
            if(this.automatic_input.cytoscape_path && this.automatic_input.cluego_base_path) {
                return "The Cytoscape executable and the ClueGO configuration files were automatically detected.  They can be changed in the section below.";
            } else if(this.automatic_input.cytoscape_path) {
                return "The Cytoscape executable was automatically detected.  It can be changed in the section below.";
            } else if(this.automatic_input.cluego_base_path) {
                return "The ClueGO configuration files were automatically detected.  It can be changed in the section below.";
            }
            return "";
        },
    },
    watch: {
        "cluego_pathways.query": function() {
            this.cluego_pathways.page = 1;
        }
    }
});