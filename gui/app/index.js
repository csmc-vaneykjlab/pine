const Vue = require("vue/dist/vue.min.js");
const { remote } = require("electron");
const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const process = require("process");
const os = require("os");
const readline = require("readline");
const shell = require("electron").shell;
const http = require("http");

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
    CREDITS: "credits",
};
const ONTOLOGY_SOURCE_TYPES = [
    {"phrase": "biologicalprocess", "name": "Biological Process"},
    {"phrase": "cellularcomponent", "name": "Cellular Component"},
    {"phrase": "molecularfunction", "name": "Molecular Function"},
    {"phrase": "human-diseases", "name": "Pathway (CLINVAR)"},
    {"phrase": "kegg", "name": "Pathway (KEGG)"},
    {"phrase": "wikipathways", "name": "Pathway (Wiki)"},
    {"phrase": "reactome", "name": "Pathway (REACTOME)"},
    {"phrase": "corum", "name": "Pathway (CORUM)"},
    {"phrase": "pathways", "name": "Pathway (Other)"},
];
const NON_NUMERIC_SORT_COLUMNS = ["GOTerm"];

const OUT_NAME_INVALID_REGEX = /[\/\\:\*\?"<>\|]/g;

function is_dir(dirname) {
    return fs.existsSync(dirname) && fs.statSync(dirname).isDirectory();
}

function is_file(filename) {
    return fs.existsSync(filename) && fs.statSync(filename).isFile();
}

const LICENSE_DIRECTORY = "app/assets/licenses";

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
        extra_required_fields: ["mods", "fasta_file", "enzyme"],
        species_map: {
            "homo sapiens": {name: "human", genemania: "4"},
            "mus musculus": {name: "mouse", genemania: "5"},
            "rattus norvegicus": {name: "rat", genemania: "7"},
            "escherichia coli": {name: "E. coli", genemania: "9"},
            "saccharomyces cerevisiae s288c": {name: "yeast", genemania: "6"},
            "arabidopsis thaliana": {name: "arabidopsis", genemania: "1"},
            "caenorhabditis elegans": {name: "C. elegans", genemania: "2"},
            "danio rerio": {name: "zebrafish", genemania: "8"},
            "drosophila melanogaster": {name: "fruit fly", genemania: "3"},
            "bos taurus": {name: "bovine", genemania: null},
            "gallus gallus": {name: "chicken", genemania: null},
            "sus scrofa": {name: "pig", genemania: null},
            "oryctolagus cuniculus": {name: "rabbit", genemania: null},
            "ovis aries": {name: "sheep", genemania: null},
            "canis lupus familiaris": {name: "dog", genemania: null},
        },
        allowed_visualize: ["biological process","cellular component","molecular function","pathways","all"],
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
            output_name: "",
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
            remove_ambiguous: null,
        },
        session_dir: null,
        missing_files: {
            in: false,
            output: false,
            fasta_file: false,
            reference_path: false,
        },
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
            labels: null,
            picked_label: null,
        },
        show_config: false,
        show_about: false,
        reanalysis_name: "",
        last_reanalysis_name: "",
        pine: null,
        credits: [
            {
                "name": "Vue.js",
                "url": "https://vuejs.org/",
                "license_file": "vue.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "Font Awesome",
                "url": "https://fontawesome.com/",
                "license_file": "font-awesome.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "Electron",
                "url": "https://electronjs.org/",
                "license_file": "electron.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "electron-builder",
                "url": "https://www.electron.build/",
                "license_file": "electron-builder.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "Sass",
                "url": "https://sass-lang.com/",
                "license_file": "sass.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "PyInstaller",
                "url": "https://www.pyinstaller.org/",
                "license_file": "pyinstaller.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "Cytoscape",
                "url": "https://cytoscape.org/",
                "license_file": "cytoscape.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "py2cytoscape",
                "url": "https://github.com/cytoscape/py2cytoscape",
                "license_file": "py2cytoscape.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "GeneMANIA",
                "url": "https://genemania.org/",
                "show_license": false,
            },
            {
                "name": "STRING",
                "url": "https://string-db.org/",
                "license_file": "string-db.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "ClueGO",
                "url": "http://www.ici.upmc.fr/cluego/",
                "license_file": "cluego.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "pandas",
                "url": "https://pandas.pydata.org/",
                "license_file": "pandas.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "cyREST",
                "url": "https://github.com/cytoscape/cyREST",
                "license_file": "cyrest.txt",
                "license": null,
                "show_license": false,
            },
            {
                "name": "igraph",
                "url": "https://igraph.org/",
                "license_file": "igraph.txt",
                "license": null,
                "show_license": false,
            },
        ],
    },
    methods: {
        run: async function(args, is_reanalysis) {
            let that = this;
            if(!this.runnable() || this.running) {
                return false;
            }
            this.running = true;

            let file_check = this.runnable_file_check(is_reanalysis);
            if(!file_check["success"]) {
                error_popup("File error", `${file_check["file"]}`);
                this.running = false;
                return false;
            }

            this.switchTab(TABS.PROGRESS);
            this.stdout = "";
            this.stderr = "";
            let stderr = "";

            this.save_settings(this.get_settings_file());

            if(process.env.NODE_ENV === "dev") {
                let args1 = [path.join(__dirname, "/../../pine/pine.py")].concat(args);
                this.pine = spawn("python", args1);
            } else {
                this.pine = spawn(path.join(__dirname, "../../extra-resources/pine/pine.exe"), args);
            }

            let new_session_dir = null;
            let cyrest_port = null;
            this.pine.stdout.on("data", function(d) {
                if(typeof d !== "string") {
                    d = d.toString("utf8");
                }
                const d_split = d.split("\n");
                for(let i = 0; i < d_split.length; i++) {
                    const ds = d_split[i].replace(/^\s+|\s+$/g, ""); // strip whitespace
                    if(ds.startsWith("COMMAND")) {
                        if(ds.startsWith("COMMAND FILE-SESSION ")) {
                            new_session_dir = ds.replace(/^COMMAND FILE-SESSION /, "");
                        }
                        if(ds.startsWith("COMMAND CYREST-PORT ")) {
                            cyrest_port = ds.replace(/^COMMAND CYREST-PORT /, "");
                        }
                    } else if(ds.length > 0 || i < d_split.length - 1) {
                        /* print every element except the last one */
                        that.stdout += ds + "\n";
                    }
                }
            });

            this.pine.stderr.on("data", function(d) {
                stderr += d + "\n";
            });

            let pr = new Promise(function(resolve, _reject) {
                that.pine.on("close", function(code) { 
                    let window = remote.getCurrentWindow();
                    if(!window.isFocused()) {
                        remote.getCurrentWindow().flashFrame(true);
                    }
                    that.running = false;
                    if(code === 0) {
                        that.stdout += "";
                        if(new_session_dir) {
                            that.new_session(new_session_dir);
                        }
                        that.switchTab(TABS.PATHWAY_SELECTION);
                        resolve(true);
                    } else {
                        that.stdout += "Run failed\n";
                        that.stderr = stderr;
                        if(cyrest_port != null) {
                            http.get(`http://localhost:${cyrest_port}/v1/commands/command/quit`);
                        }

                        if(!is_reanalysis) {
                            /* delete the orphaned directory if it exists */
                            if(new_session_dir && is_dir(new_session_dir)) {
                                const dir_contents = fs.readdirSync(new_session_dir);
                                if(dir_contents.length == 2 && dir_contents.includes("PINE.log") && dir_contents.includes("timestamp.json")) {
                                    const pine_log_file = path.join(new_session_dir, "PINE.log");
                                    const pine_settings_file = path.join(new_session_dir, "timestamp.json");
                                    if(is_file(pine_log_file) && is_file(pine_settings_file)) {
                                        fs.unlinkSync(pine_log_file);
                                        fs.unlinkSync(pine_settings_file);
                                        fs.rmdirSync(new_session_dir);
                                    }
                                }
                            }
                        }
                        resolve(false);
                        document.getElementById("log-wrapper").scrollTop = 0;
                    }
                });
            });

            return await pr;
        },
        run_full: function() {
            /* if genemania can't be found, show popup confirming they want to continue */
            const gm_check = this.genemania_check();
            if(gm_check != null) {
                const res = remote.dialog.showMessageBoxSync({
                    "type": "warning",
                    "buttons": [
                        "Continue anyways",
                        "Cancel (recommended)",
                    ],
                    "defaultId": 1,
                    "title": "GeneMANIA error",
                    "message": gm_check,
                });
                if(res !== 0) {
                    return;
                }
            }

            let args = this.pine_args();
            this.run(args, false);
        },
        run_with_cluego_subset: async function() {
            this.input.output_name = "";
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

            if(!is_dir(this.session_dir)) {
                error_popup("Invalid session directory", "Your session directory may have been renamed or deleted.  Please reload the session or run this dataset again.", true);
                return;
            }

            if(!this.unique_reanalysis_name(reanalysis_name)) {
                error_popup("Name in use", "This name has already been used for this session, please pick a different name", true);
                return;
            }
            this.last_reanalysis_name = reanalysis_name;
            const filtered_file_name = path.join(this.session_dir, reanalysis_name + ".cluego.txt");
            const log_file_name = path.join(this.session_dir, reanalysis_name + ".log");
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
            let res = await this.run(args, true);
            if(!res) {
                if(is_file(filtered_file_name)) {
                    fs.unlinkSync(filtered_file_name);
                }
                if(is_file(log_file_name)) {
                    fs.unlinkSync(log_file_name);
                }
            }
        },
        genemania_check: function() {
            if(this.input.run !== "both" && this.input.run !== "genemania") {
                return null; // check passed because genemania is not needed
            }

            const msg_genemania = "Please install Genemania plugin. If it is already installed, ignore this message and continue.";
            const msg_invalid_species = "Unsupported species chosen";
            const msg_missing_species = `Please install '${this.input.species}' dataset for Genemania plugin within Cytoscape. If it is already installed, ignore this message and continue.`;

            /* check for genemania configuration directory - only should need to be checked on first run */
            let gm_config_dir = path.join(os.homedir(), "Documents/genemania_plugin");
            if(!is_dir(gm_config_dir)) {
                return msg_genemania;
            }

            const files = fs.readdirSync(gm_config_dir);
            let subdirs = [];
            for(const f of files) {
                const subdir_full_path = path.join(gm_config_dir, f);
                if(!is_dir(subdir_full_path)) {
                    continue;
                }
                if(f.startsWith("gmdata-")) {
                    subdirs.push(subdir_full_path);
                }
            }
            if(subdirs.length === 0) {
                return msg_genemania;
            }

            const species_number = this.get_genemania_species(this.input.species);
            if(species_number == null) {
                return msg_invalid_species;
            }

            for(const sd of subdirs) {
                const sd_files = fs.readdirSync(sd);
                for(const sdf of sd_files) {
                    if(!is_dir(path.join(sd, sdf))) {
                        continue;
                    }
                    if(sdf === species_number && is_file(path.join(sd, sdf, "metadata.xml"))) {
                        return null;
                    }
                }
            }

            return msg_missing_species;
        },
        cancel_pine: function() {
            if(this.pine === null || !this.running) {
                return;
            }
            this.pine.kill();
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
                pad(now.getMonth() + 1) +
                pad(now.getDate()) +
                "_" +
                pad(now.getHours()) +
                pad(now.getMinutes()) +
                pad(now.getSeconds());
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
        },
        reset_session: function() {
            this.session_dir = null;
            this.input.output_name = "";
            this.reset_cluego_pathways();
            this.switchTab(TABS.INPUT);
        },
        set_input: function(name, val) {
            this.input[name] = val;
        },
        session_exists: function() {
            if(this.session_dir) {
                return true;
            }
            return false;
        },
        user_load_session: function(e) {
            if(!e.target.files[0].path) {
                e.target.value = "";
                return;
            }
            let new_dir = e.target.files[0].path;
            if(!is_dir(new_dir)) {
                error_popup("Invalid path", "Path provided is not a directory");
                e.target.value = "";
                return;
            }
            let old_dir = this.session_dir;
            this.session_dir = new_dir;
            if(!is_file(this.session_cluego_file) || !is_file(this.session_settings_file) || !is_file(this.session_timestamp_file)) {
                this.session_dir = old_dir;
                e.target.value = "";
                error_popup("Invalid session", "The session directory you provided is not valid");
                return;
            }
            this.load_settings(this.session_settings_file);
            const file_check = this.runnable_file_check(true);
            if(!file_check["success"]) {
                error_popup("Invalid session", `Input file ${file_check["filename"]} has been moved or deleted.  Please restore this file to original location or start a new session.`);
            }
            this.read_cluego_pathways();
            if(file_check["success"]) {
                this.switchTab(TABS.PATHWAY_SELECTION);
            }
            e.target.value = "";
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
                "--cytoscape-exe", this.input.cytoscape_path,
                "--gui",
            ];
            if(this.is_extra_options_required()) {
                args.push("--mods");
                args.push(this.input.mods);

                args.push("--enzyme");
                args.push(this.input.enzyme);

                args.push("--fasta-file");
                args.push(this.input.fasta_file);
            }
            if(this.input.output_name) {
                args.push("--output-name");
                args.push(this.input.output_name);
            }
            if(this.input.significant) {
                args.push("--significant");
            }
            if(this.input.remove_ambiguous) {
                args.push("--exclude-ambiguity");
            }

            return args;
        },
        runnable: function() {
            const input_check =
                this.input.in &&
                this.get_cluego_mapping() &&
                this.input.output &&
                this.input.type &&
                !this.running &&
                this.validate_inputs();
            if(input_check) {
                if(this.is_extra_options_required()) {
                    for(const erf of this.extra_required_fields) {
                        if(!this.input[erf]) {
                            return false;
                        }
                    }
                    return true;
                } else {
                    return true;
                }
            }
            return false;
        },
        runnable_file_check: function(is_reanalysis) {
            let msg = null;
            let filename = null;
            this.reset_missing_files();
            if(!this.input.in || this.missing_file_check('in', false)) {
                msg = "Input file doesn't exist.";
                filename = this.input.in;
            } else if(!this.input.output || this.missing_file_check('output', true)) {
                msg = "Output directory doesn't exist.";
                filename = this.input.output;
            } else if(this.is_extra_options_required() && (!this.input.fasta_file || this.missing_file_check('fasta_file', false))) {
                msg = "Fasta file doesn't exist.";
                filename = this.input.fasta_file;
            } else if(this.input.reference_path && this.missing_file_check('reference_path', false)) {
                msg = "Reference path file doesn't exist.";
                filename = this.input.reference_path;
            } else if(this.input.output_name && !is_reanalysis) {
                const output_name_path = path.join(this.input.output, this.input.output_name);
                if(is_file(output_name_path) || is_dir(output_name_path)) {
                    msg = `Output path ${output_name_path} already exists.`;
                }
                filename = output_name_path;
            }
            if(msg) {
                return {
                    "success": false,
                    "file": msg,
                    "filename": filename,
                };
            }
            return {"success": true};
        },
        missing_file_check: function(var_name, check_is_dir) {
            if(!this.input[var_name]) {
                return this.missing_files[var_name] = false;
            }
            let exists;
            if(check_is_dir) {
                exists = is_dir(this.input[var_name]);
            } else {
                exists = is_file(this.input[var_name]);
            }
            return this.missing_files[var_name] = !exists;
        },
        reset_missing_files: function() {
            for(var_name in this.missing_files) {
                this.missing_files[var_name] = false;
            }
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
            this.cluego_pathways.labels = [];
            this.cluego_pathways.picked_label = null;
        },
        validate_inputs: function() {
            if(this.is_extra_options_required() && !this.validate_inputs_mods()) {
                return false;
            } else if(
                !this.validate_inputs_in() ||
                !this.validate_inputs_fasta_file() ||
                !this.validate_inputs_reference_path() ||
                !this.validate_inputs_fccutoff() ||
                !this.validate_inputs_pvalcutoff() ||
                !this.validate_inputs_score() ||
                !this.validate_inputs_limit() ||
                !this.validate_inputs_cluego_pval()
            ) {
                return false;
            }
            return true;
        },
        validate_inputs_in: function() {
            if(!this.input.in) {
                return true;
            }
            return this.input.in.endsWith(".csv");
        },
        validate_inputs_fasta_file: function() {
            if(!this.input.fasta_file) {
                return true;
            }
            return this.input.fasta_file.endsWith(".fasta");
        },
        validate_inputs_reference_path: function() {
            if(!this.input.reference_path) {
                return true;
            }
            return this.input.reference_path.endsWith(".txt");
        },
        validate_inputs_mods: function() {
            /* allow for a comma separated list of amino acids with modifications */
            const allowed_chars = "[^[\\](){}]+";
            const single_element = `[A-Z](?:\\(${allowed_chars}\\)|\\[${allowed_chars}\\]|\\{${allowed_chars}\\})?`;
            const regex = RegExp(`^${single_element}(?:, ?${single_element})*$`)
            return regex.test(this.input.mods);
        },
        validate_inputs_fccutoff: function() {
            const parsed = parseFloat(this.input.fccutoff);
            if(isNaN(parsed)) {
                return false;
            }
            return parsed >= 0.0;
        },
        validate_inputs_pvalcutoff: function() {
            const parsed = parseFloat(this.input.pvalcutoff);
            if(isNaN(parsed)) {
                return false;
            }
            return parsed >= 0.0 && parsed <= 1.0;
        },
        validate_inputs_score: function() {
            const parsed = parseFloat(this.input.score);
            if(isNaN(parsed)) {
                return false;
            }
            return parsed >= 0.0 && parsed <= 0.99;
        },
        validate_inputs_limit: function() {
            const parsed = parseInt(this.input.limit);
            if(isNaN(parsed)) {
                return false;
            }
            return parsed >= 0 && parsed <= 100;
        },
        validate_inputs_cluego_pval: function() {
            const parsed = parseFloat(this.input.cluego_pval);
            if(isNaN(parsed)) {
                return false;
            }
            return parsed >= 0.0 && parsed <= 1.0;
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
            this.input.pvalcutoff = 0.05;
            this.input.visualize = "pathways";
            this.input.cluego_pval = 0.05;
            this.input.reference_path = "";
            this.reset_input_file("reference_path");
            this.input.grouping = "medium";
            this.input.enzyme = "trypsin";
            this.reset_input_file("fasta_file");
            this.input.mods = "S,T,Y";
            this.reset_input_file("in");
            this.reset_input_file("output");
            this.input.output_name = "";
            this.input.remove_ambiguous = true;
        },
        reset_input_file: function(name) {
            const refs_name = "input_" + name;

            this.input[name] = "";
            if(this.$refs.hasOwnProperty(refs_name)) {
                this.$refs[refs_name].value = "";
            }
        },
        read_cluego_pathways: function() {
            var that = this;

            this.reset_cluego_pathways();

            if(!this.session_exists() || !is_file(this.session_cluego_file)) {
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
                    let record = {"data": {}, "selected": false, "line": line, "labels": {}};
                    for(let i = 0; i < that.cluego_pathways.header.length; i++) {
                        const header = that.cluego_pathways.header[i];
                        if(header.startsWith("% change")) {
                            let label;
                            if(header.includes(":")) {
                                label = header.split(":").slice(1).join(":");
                            } else {
                                label = "Status";
                            }
                            const val = parseFloat(fields[i]);
                            if(isNaN(val)) {
                                record["labels"][label] = 0.0;
                            } else {
                                record["labels"][label] = val;
                            }
                        }
                        record["data"][header] = fields[i];
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

            line_reader.on("close", () => {
                this.cluego_pathways.data = cluego_pathways;

                /* get labels */
                let labels = new Set();
                for(const record of this.cluego_pathways.data) {
                    for(const label in record.labels) {
                        labels.add(label);
                    }
                }
                this.cluego_pathways.labels = Array.from(labels);
                this.cluego_pathways.labels.sort((x, y) => x.localeCompare(y));
                this.cluego_pathways.picked_label = this.cluego_pathways.labels.length > 0 ? this.cluego_pathways.labels[0] : "";
            });
        },
        get_cluego_mapping: function() {
            if(!this.cluego_picked_version) {
                return null;
            }

            for(const mapping of this.cluego_picked_version.mapping_files) {
                if(this.species_map[mapping.species.toLowerCase()].name === this.input.species) {
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
                this.setCytoscapePath(path, true);
            } else if(name === "cluego_base_path") {
                this.setCluegoBasePath(path, true);
            }

            this.input[name] = path;

            this.runnable_file_check(false);
        },
        save_settings: function(settings_file) {
            let settings = {};
            for(const inp in this.input) { 
                if(inp === "output_name") {
                    continue;
                }
                if(!this.is_extra_options_required() && this.extra_required_fields.includes(inp)) {
                    continue;
                }
                settings[inp] = this.input[inp];
            }
            settings.cluego_picked_version = this.cluego_picked_version.version;
            fs.writeFileSync(settings_file, JSON.stringify(settings, null, 4));
        },
        load_settings: function(settings_file) {
            if(!is_file(settings_file)) {
                return;
            }
            const raw = fs.readFileSync(settings_file, "utf8");
            let cluego_picked_version = null;
            try {
                const saved_input = JSON.parse(raw);
                for(const key in saved_input) {
                    if(key in this.input) {
                        if(key === "cluego_base_path") {
                            this.setCluegoBasePath(saved_input[key], false);
                        } else if(key === "cytoscape_path") {
                            this.setCytoscapePath(saved_input[key], false);
                        } else {
                            Vue.set(this.input, key, saved_input[key]);
                        }
                    } else if(key === "cluego_picked_version") {
                        cluego_picked_version = saved_input[key];
                    }
                }
            } catch(e) {
            }

            if(cluego_picked_version != null) {
                for(const cv of this.cluego_versions) {
                    if(cv.version === cluego_picked_version) {
                        this.cluego_picked_version = cv;
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
                            this.setCytoscapePath(cyto_path, false);
                            break;
                        }
                    }
                }
            }

            if(!this.input.cluego_base_path) {
                const cluego_base_path = path.join(os.homedir(), "ClueGOConfiguration");
                if(fs.existsSync(cluego_base_path)) {
                    this.setCluegoBasePath(cluego_base_path, false);
                }
            }

        },
        setCluegoBasePath: function(cluego_path, show_warnings) {
            if(!fs.existsSync(cluego_path)) {
                if(show_warnings) {
                    error_popup("Path does not exist", "ClueGO configuration file path does not exist");
                }
                return;
            }

            const path_split = cluego_path.split(path.sep);
            const path_index = path_split.indexOf(CLUEGO_CONFIGURATION_BASE_NAME);
            if(path_index === -1) {
                if(show_warnings) {
                    error_popup("Invalid directory", "ClueGO configuration directory must be named " + CLUEGO_CONFIGURATION_BASE_NAME);
                }
                return;
            }

            cluego_path = path_split.slice(0, path_index + 1).join(path.sep);
            if(!fs.statSync(cluego_path).isDirectory()) {
                if(show_warnings) {
                    error_popup("Invalid directory", CLUEGO_CONFIGURATION_BASE_NAME + " must be a directory");
                }
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
                if(show_warnings) {
                    error_popup("Invalid configuration directory", "This directory did not any valid accession files.");
                }
            }
        },
        setCytoscapePath: function(cy_path, show_warnings) {
            if(!is_file(cy_path)) {
                if(show_warnings) {
                    error_popup("Invalid path", "Cytoscape.exe must be a file");
                }
                return;
            }
            if(!cy_path.toLowerCase().endsWith("cytoscape.exe")) {
                if(show_warnings) {
                    error_popup("Invalid file", "Please provide the Cytoscape.exe executable file");
                }
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
                    return this.running || this.stdout;
                case TABS.PATHWAY_SELECTION:
                    return !this.running && this.session_exists();
                case TABS.CREDITS:
                    return true;
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
        open_url: function(url) {
            shell.openExternal(url);
        },
        toggle_pathway_selected: function(pathway) {
            pathway.selected = !pathway.selected;
        },
        toggle_show_license: function(credit) {
            if(credit.show_license) {
                credit.show_license = false;
                return;
            }

            if(credit.license == null) { // double equals catches null and undefined
                if(credit.license_file == null) {
                    return;
                }
                fs.readFile(path.join(LICENSE_DIRECTORY, credit.license_file), "utf-8", (err, data) => {
                    if(err) {
                        console.log(err);
                        return;
                    }
                    credit.license = data;
                    credit.show_license = true;
                });
            } else {
                credit.show_license = true;
            }
        },
        pathway_data_label: function(datum, label) {
            let percent = datum.labels[label];

            if(percent == null || isNaN(percent)) {
                return NaN;
            }

            let icon;
            if(percent > 0) {
                icon = {
                    "classes": "fas fa-arrow-up color-up-reg",
                    "tooltip": "Upregulation",
                };
            } else if(percent < 0) {
                icon = {
                    "classes": "fas fa-arrow-down color-down-reg",
                    "tooltip": "Downregulation",
                };
            } else {
                icon = {
                    "classes": "fas fa-minus",
                    "tooltip": "No change",
                };
            }
            let display;
            if(percent > 0) {
                display = `${percent}%`;
            } else if(percent < 0) {
                display = `${-percent}%`;
            } else {
                display = "N/A";
            }

            return `
                <span class="tooltip-parent">
                    <i class="${icon.classes}"></i>
                    <div class="tooltip">${icon.tooltip}</div>
                </span>
                ${display}
            `;
        },
        get_genemania_species: function(species_name) {
            for(const species_key in this.species_map) {
                const species = this.species_map[species_key];
                if(species.name === species_name) {
                    return species.genemania;
                }
            }
            return null;
        },
        selected_species_updated: function(species_name) {
            if(this.get_genemania_species(species_name) == null) {
                this.input.run = "string";
            }
        },
    },
    mounted: function() {
        this.reset_cluego_pathways();
        this.set_input_defaults();
        this.load_settings(this.get_settings_file());
        this.searchForPaths();
        this.credits.sort((a, b) => {
            if(a.name.toLowerCase() < b.name.toLowerCase()) {
                return -1;
            }
            if(a.name.toLowerCase() > b.name.toLowerCase()) {
                return 1;
            }
            return 0;
        });
    },
    filters: {
        longname: function(v, len, is_path) {
            if(v == null) {
                return "";
            }
            let basename;
            if(is_path) {
                basename = path.basename(v);
            } else {
                basename = v;
            }
            if(basename.length > len) {
                return basename.slice(0, len - 3) + "...";
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
        },
        ontology_source_category_name: function(record) {
            if(record["Ontology Source Category"] === "Other") {
                return record["Ontology Source"];
            }
            return record["Ontology Source Category"];
        },
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
                } else if(col === "label") {
                    filtered.sort(function(a, b) {
                        if(a.labels[that.cluego_pathways.picked_label] < b.labels[that.cluego_pathways.picked_label]) {
                            return lower;
                        } else if(a.labels[that.cluego_pathways.picked_label] > b.labels[that.cluego_pathways.picked_label]) {
                            return higher;
                        } else {
                            return 0;
                        }
                    });
                } else {
                    filtered.sort(function(a, b) {
                        if(!(col in a.data)) return 1;
                        if(!(col in b.data)) return 0;
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
                labels: this.cluego_pathways.labels,
                picked_label: this.cluego_pathways.picked_label,
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
            return path.join(this.session_dir, "settings-gui.json");
        },
        session_timestamp_file: function() {
            if(!this.session_dir) {
                return "";
            }
            return path.join(this.session_dir, "timestamp.json");
        },
        selectable_species: function() {
            let selectable = [];
            let seen = new Set();
            for(const mapping of this.cluego_picked_version.mapping_files) {
                const mapping_species = mapping.species.toLowerCase();
                if(mapping_species in this.species_map) {
                    const name = this.species_map[mapping_species].name;
                    let display_name = name;
                    if(this.species_map[mapping_species].genemania == null) {
                        display_name += " (STRING only)";
                    }
                    selectable.push({
                        "name": name,
                        "display_name": display_name,
                        "selectable": true,
                        "has_genemania": this.species_map[mapping_species].genemania != null,
                    });
                    seen.add(this.species_map[mapping_species].name);
                }
            }
            for(const key in this.species_map) {
                const species = this.species_map[key].name;
                if(seen.has(species)) {
                    continue;
                }
                selectable.push({
                    "name": species,
                    "display_name": species + " (not installed in ClueGO)",
                    "selectable": false,
                });
                seen.add(species);
            }
            return selectable;
        },
        allowed_runs: function() {
            let runs = [
                {"name": "string", "disabled": false},
                {"name": "genemania", "disabled": false},
                {"name": "both", "disabled": false},
            ];
            if(this.get_genemania_species(this.input.species) == null) {
                runs = runs.map((x) => {
                    if(x.name === "string") {
                        return x;
                    } else {
                        return Object.assign({}, x, {"disabled": true});
                    }
                });
            }
            return runs;
        },
    },
    watch: {
        "cluego_pathways.query": function() {
            this.cluego_pathways.page = 1;
        },
        "cluego_pathways.ontology_sources_filter": function() {
            this.cluego_pathways.page = 1;
        },
        "input.output_name": function(new_val) {
            this.$set(this.input, "output_name", new_val.replace(OUT_NAME_INVALID_REGEX, ""));
        },
        "reanalysis_name": function(new_val) {
            this.$set(this, "reanalysis_name", new_val.replace(OUT_NAME_INVALID_REGEX, ""));
        },
    }
});