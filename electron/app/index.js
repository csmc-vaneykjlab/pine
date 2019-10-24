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

const CLUEGO_CONFIGURATION_BASE_NAME = "ClueGOConfiguration";
const TABS = {
    SETUP: "setup",
    INPUT: "input",
    PROGRESS: "progress",
    PATHWAY_SELECTION: "pathway-selection",
};
const ONTOLOGY_SOURCE_TYPES = [
    {"phrase": "biologicalprocess", "name": "Biological Process"},
    {"phrase": "cellularcomponent", "name": "Cellular Component"},
    {"phrase": "molecularfunction", "name": "Molecular Function"},
    {"phrase": "human-diseases", "name": "Pathway (Human Diseases)"},
    {"phrase": "kegg", "name": "Pathway (KEGG)"},
    {"phrase": "wikipathways", "name": "Pathway (Wiki Pathways)"},
    {"phrase": "pathways", "name": "Pathway (Other)"},
];
const NON_NUMERIC_SORT_COLUMNS = ["GOTerm"];

function is_dir(dirname) {
    return fs.existsSync(dirname) && fs.statSync(dirname).isDirectory();
}

function is_file(filename) {
    return fs.existsSync(filename) && fs.statSync(filename).isFile();
}

function error_popup(title, message, warning) {
    let type = "error";
    if(warning) {
        type = "warning";
    }
    remote.dialog.showMessageBox({"type": type, "buttons": ["Ok"], "defaultId": 0, "title": title, "message": message});
}

let vm = new Vue({
    el: "#app",
    data: {
        allowed_types: {
            "noFC": {"text": "No fold change", "er": false}, // er - extra required - extra fields need to be shown and submitted by the user
            "singleFC": {"text": "Single fold change", "er": false},
            "multiFC": {"text": "Multi fold change", "er": false},
            "category": {"text": "Category", "er": false},
            "singleFC-ptm": {"text": "Single fold change PTM", "er": true},
            "multiFC-ptm": {"text": "Multi fold change PTM", "er": true},
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
        session_dir: null,
        stdout: "",
        stderr: "",
        running: false,
        tabs: TABS,
        current_tab: TABS.SETUP,
        cluego_versions: null,
        cluego_picked_version: null,
        cluego_pathways: {
            data: null,
            header: null,
            query: null,
            page: null,
            per_page: null,
            sort: null,
            ontology_sources_filter: null,
        },
        show_config: false,
        show_about: false,
        reanalysis_name: "",
        last_reanalysis_name: "",
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

            this.save_settings(this.get_settings_file());

            if(process.env.NODE_ENV === "dev") {
                let args1 = [path.join(__dirname, "/../../changes_to_pine_final.py")].concat(args);
                var pine = spawn("C:/Users/GoJ1/AppData/Local/Programs/Python/Python37/python.exe", args1);
            } else {
                var pine = spawn(path.join(__dirname, "/../../pine_2/pine_2.exe"), args);
            }

            let new_session_dir = null;
            pine.stdout.on("data", function(d) {
                if(typeof d !== "string") {
                    d = d.toString("utf8");
                }
                const d_split = d.split("\n");
                for(let ds of d_split) {
                    ds = ds.replace(/^\s+|\s+$/g, ""); // strip whitespace
                    if(ds.startsWith("COMMAND")) {
                        if(ds.startsWith("COMMAND FILE-SESSION ")) {
                            new_session_dir = ds.replace(/^COMMAND FILE-SESSION /, "");
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
                    if(new_session_dir) {
                        that.new_session(new_session_dir);
                    }
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
            if(this.reanalysis_name) {
                var reanalysis_name = this.reanalysis_name;
                this.reanalysis_name = "";
            } else {
                var reanalysis_name = this.generate_reanalysis_name();
            }
            if(!this.unique_reanalysis_name(reanalysis_name)) {
                error_popup("Name in use", "This name has already been used for this session, please pick a different name", true);
                return;
            }
            this.last_reanalysis_name = reanalysis_name;
            const filtered_file_name = path.join(this.session_dir, reanalysis_name + ".cluego.txt");
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
            this.run(args);
        },
        generate_reanalysis_name: function() {
            function pad(n) {
                if(n < 10) {
                    return "0" + n.toString();
                }
                return n.toString();
            }
            let now = new Date();
            let timestamp =
                now.getFullYear().toString() +
                pad(now.getMonth()) +
                pad(now.getDate()) +
                "T" +
                pad(now.getHours()) +
                pad(now.getMinutes());
            return timestamp + "_reanalysis_PINE";
        },
        unique_reanalysis_name: function(name) {
            if(!is_dir(this.session_dir)) {
                return false;
            }
            check_files = [
                path.join(this.session_dir, name + ".cluego.txt"),
                path.join(this.session_dir, name + ".cys"),
                path.join(this.session_dir, name + ".log"),
            ];
            for(const cf of check_files) {
                if(fs.existsSync(cf)) {
                    return false;
                }
            }
            return true;
        },
        new_session: function(dir) {
            if(!dir || !is_dir(dir)) {
                error_popup("Session not created", "The new session could not be created");
                return;
            }
            this.session_dir = dir;
            if(!this.session_cluego_file || !this.session_cytoscape_file || !this.session_settings_file) {
                this.session_dir = null;
            }
            this.save_settings(this.session_settings_file);
            this.read_cluego_pathways();
            this.switchTab(TABS.PATHWAY_SELECTION);
        },
        reset_session: function() {
            this.session_dir = null;
            this.reset_cluego_pathways();
            this.switchTab(TABS.INPUT);
        },
        session_exists: function() {
            if(this.session_dir) {
                return true;
            }
            return false;
        },
        user_load_session: function(e) {
            if(!e.target.files[0].path) {
                return;
            }
            let new_dir = e.target.files[0].path;
            if(!is_dir(new_dir)) {
                error_popup("Invalid path", "Path provided is not a directory");
                return;
            }
            let old_dir = this.session_dir;
            this.session_dir = new_dir;
            if(!is_file(this.session_cluego_file) || !is_file(this.session_cytoscape_file) || !is_file(this.session_settings_file)) {
                this.session_dir = old_dir;
                error_popup("Invalid session", "The session directory you provided is not valid");
            }
            this.load_settings(this.session_settings_file);
            this.read_cluego_pathways();
            this.switchTab(TABS.PATHWAY_SELECTION);
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
            this.cluego_pathways.sort = null;
            this.cluego_pathways.ontology_sources_filter = "All";
        },
        set_input_defaults: function() {
            if(!this.settings_editable()) {
                return;
            }
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

            if(!this.session_exists() || !fs.existsSync(this.session_cluego_file) || !fs.statSync(this.session_cluego_file).isFile()) {
                return;
            }

            let counter = 0;
            let cluego_pathways = [];
            let line_reader =  readline.createInterface({
                input: fs.createReadStream(this.session_cluego_file),
            });
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

                    /* get ontology source category */
                    record["data"]["Ontology Source Category"] = "Other";   // initialize to other
                    for(const st of ONTOLOGY_SOURCE_TYPES) {
                        if(record["data"]["Ontology Source"].toLowerCase().includes(st.phrase)) {
                            record["data"]["Ontology Source Category"] = st.name;
                            break;  // record should only match a single ontology source
                        }
                    }

                    cluego_pathways.push(record);
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
            } else if(name === "cluego_base_path") {
                this.setCluegoBasePath(path);
            } else {
                this.input[name] = path;
            }
            this.refreshTab();
        },
        save_settings: function(settings_file) {
            fs.writeFileSync(settings_file, JSON.stringify(this.input, null, 4));
        },
        load_settings: function(settings_file) {
            if(!is_file(settings_file)) {
                return;
            }
            const raw = fs.readFileSync(settings_file, "utf8");
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
        get_settings_file: function() {
            let app_dir = process.env.APPDATA || (process.platform == 'darwin' ? process.env.HOME + 'Library/Preferences' : process.env.HOME + "/.local/share");
            let full_dir = path.join(app_dir, "Cedars Sinai JVE", "PINE");
            if(!fs.existsSync(full_dir)) {
                fs.mkdirSync(full_dir, {recursive: true});
            }
            return path.join(full_dir, "settings.json");
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
                            break;
                        }
                    }
                }
            }

            if(!this.input.cluego_base_path) {
                const cluego_base_path = path.join(os.homedir(), "ClueGOConfiguration");
                if(fs.existsSync(cluego_base_path)) {
                    this.setCluegoBasePath(cluego_base_path);
                }
            }

            this.refreshTab();
        },
        refreshTab: function() {
            if(this.current_tab === TABS.SETUP && this.input.cytoscape_path && this.input.cluego_base_path) {
                this.switchTab(TABS.INPUT);
            }
        },
        setCluegoBasePath: function(cluego_path) {
            if(!fs.existsSync(cluego_path)) {
                error_popup("Path does not exist", "ClueGO configuration file path does not exist");
                return;
            }

            const path_split = cluego_path.split(path.sep);
            const path_index = path_split.indexOf(CLUEGO_CONFIGURATION_BASE_NAME);
            if(path_index === -1) {
                error_popup("Invalid directory", "ClueGO configuration directory must be named " + CLUEGO_CONFIGURATION_BASE_NAME);
                return;
            }

            cluego_path = path_split.slice(0, path_index + 1).join(path.sep);
            if(!fs.statSync(cluego_path).isDirectory()) {
                error_popup("Invalid directory", CLUEGO_CONFIGURATION_BASE_NAME + " must be a directory");
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
            } else {
                error_popup("Invalid configuration directory", "This directory did not any valid accession files.");
            }
        },
        setCytoscapePath: function(cy_path) {
            if(!is_file(cy_path)) {
                error_popup("Invalid path", "Cytoscape.exe must be a file");
                return;
            }
            if(!cy_path.toLowerCase().endsWith("cytoscape.exe")) {
                error_popup("Invalid file", "Please provide the Cytoscape.exe executable file");
                return;
            }
            this.input.cytoscape_path = cy_path;
        },
        open_last_cytoscape_file: function() {
            if(!this.session_dir) {
                return;
            }
            if(this.last_reanalysis_name) {
                var cys_file = path.join(this.session_dir, this.last_reanalysis_name + ".cys");
            } else {
                var cys_file = path.join(this.session_dir, "PINE.cys");
            }
            if(!is_file(cys_file)) {
                error_popup("File doesn't exist", "Could not find last Cytoscape file for this session.");
                return;
            }
            shell.openItem(cys_file);
        },
        open_session_dir: function() {
            if(!this.session_dir || !is_dir(this.session_dir)) {
                error_popup("Path doesn't exist", "Session is not set or session directory does not exist.");
                return;
            }
            shell.openItem(this.session_dir);
        },
        tab_reachable: function(tab_name) {
            switch(tab_name) {
                case TABS.SETUP:
                    return true;
                case TABS.INPUT:
                    return this.input.cytoscape_path && this.input.cluego_base_path;
                case TABS.PROGRESS:
                    return this.running || this.stdout || this.session_exists();
                case TABS.PATHWAY_SELECTION:
                    return !this.running && this.session_exists();
            }
            return false;
        },
        settings_editable: function() {
            return !this.running && !this.session_exists();
        },
        switchTab: function(tab_name) {
            if(this.tab_reachable(tab_name)) {
                this.current_tab = tab_name;
            }
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
        cluego_pathways_sort: function(column, display) {
            let numeric = NON_NUMERIC_SORT_COLUMNS.includes(column) ? false : true;
            let sort_dir = "asc";
            if(this.cluego_pathways.sort && this.cluego_pathways.sort.column === column && this.cluego_pathways.sort.order === "asc") {
                sort_dir = "desc";
            }
            this.cluego_pathways.sort = {"column": column, "order": sort_dir, "numeric": numeric, "display": display};
        },
        cluego_pathways_sort_clear: function() {
            this.cluego_pathways.sort = null;
        },
        cluego_pathways_unselect_all: function() {
            for(const cpd of this.cluego_pathways.data) {
                if(cpd.selected) {
                    cpd.selected = false;
                }
            }
        },
        validate_reanalysis_name: function() {
            let match = this.reanalysis_name.match(/[A-Za-z0-9\-_ ]/g);
            if(match) {
                this.reanalysis_name = match.join("");
            } else {
                this.reanalysis_name = "";
            }
        },
        open_url: function(url) {
            shell.openExternal(url);
        },
    },
    mounted: function() {
        this.reset_cluego_pathways();
        this.set_input_defaults();
        this.load_settings(this.get_settings_file());
        this.searchForPaths();
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

            /* filter on query */
            let filtered = this.cluego_pathways.data.filter(function(pathway) {
                if(that.cluego_pathways.query !== "") {
                    let lower = that.cluego_pathways.query.toLowerCase();
                    if(!pathway.data['GOTerm'].toLowerCase().includes(lower)) {
                        return false;
                    }
                }

                return true;
            });

            /* get ontology sources counts */
            let pathways_source_name = "Pathway (All)";
            let sources = {"All": filtered.length};
            for(const f of filtered) {
                const source_name = f.data["Ontology Source Category"];
                if(!(source_name in sources)) {
                    sources[source_name] = 0;
                }
                sources[source_name] += 1;

                if(source_name.toLowerCase().includes("pathway")) {
                    if(!(pathways_source_name in sources)) {
                        sources[pathways_source_name] = 0
                    }
                    sources[pathways_source_name] += 1
                }
            }
            let ontology_sources = [];
            for(const s in sources) {
                ontology_sources.push({"name": s, "count": sources[s]});
            }
            ontology_sources.sort(function(a,b) {
                if(a.name === "All") return -1;
                if(b.name === "All") return 1;
                if(a.name < b.name) return -1;
                if(a.name > b.name) return 1;
                return 0;
            });

            /* filter on ontology source - filter sources after getting the category counts */
            filtered = filtered.filter(function(pathway) {
                if(that.cluego_pathways.ontology_sources_filter === "All") {
                    return true;
                }

                if(that.cluego_pathways.ontology_sources_filter === pathways_source_name) {
                    return pathway.data["Ontology Source Category"].toLowerCase().includes("pathway");
                }

                return pathway.data["Ontology Source Category"] === that.cluego_pathways.ontology_sources_filter;
            });

            if(this.cluego_pathways.sort) {
                const lower = this.cluego_pathways.sort.order === "asc" ? -1 : 1;
                const higher = lower * -1;
                const col = this.cluego_pathways.sort.column;
                const numeric = this.cluego_pathways.sort.numeric;
                if(col === "selected") {
                    filtered.sort(function(a, b) {
                        if(a.selected && !b.selected) return lower;
                        if(!a.selected && b.selected) return higher;
                        return 0;
                    });
                } else {
                    filtered.sort(function(a, b) {
                        if(!(col in a.data)) return 1;
                        if(!(col in b.data)) return -1;
                        let a1 = a.data[col];
                        let b1 = b.data[col];
                        if(numeric) {
                            a1 = Number(a1);
                            if(isNaN(a1)) return 1;
                            b1 = Number(b1);
                            if(isNaN(b1)) return -1;
                        }
                        if(a1 < b1) return lower;
                        if(a1 > b1) return higher;
                        return 0;
                    });
                }
            }

            /* get the number of selected records */
            let n_selected = 0;
            for(const cpd of this.cluego_pathways.data) {
                if(cpd.selected) {
                    n_selected += 1;
                }
            }

            const full_length = filtered.length;
            const n_pages = Math.ceil(full_length / this.cluego_pathways.per_page);

            const start_idx = (this.cluego_pathways.page - 1) * this.cluego_pathways.per_page;
            const end_idx = start_idx + this.cluego_pathways.per_page;

            const data = filtered.slice(start_idx, end_idx);

            return {
                data: data,
                n_empty_rows: that.cluego_pathways.per_page - data.length,
                full_length: full_length,
                n_pages: n_pages,
                ontology_sources: ontology_sources,
                n_selected: n_selected,
            };
        },
        session_cluego_file: function() {
            if(!this.session_dir) {
                return "";
            }
            return path.join(this.session_dir, "PINE.cluego.txt");
        },
        session_cytoscape_file: function() {
            if(!this.session_dir) {
                return "";
            }
            return path.join(this.session_dir, "PINE.cys");
        },
        session_settings_file: function() {
            if(!this.session_dir) {
                return "";
            }
            return path.join(this.session_dir, "settings.txt");
        },
    },
    watch: {
        "cluego_pathways.query": function() {
            this.cluego_pathways.page = 1;
        },
        "cluego_pathways.ontology_sources_filter": function() {
            this.cluego_pathways.page = 1;
        },
    }
});