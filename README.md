# PINE Installation and Usage

## Table of Contents
- [Requirements and Setup](#requirements-and-setup)
- [Using PINE](#using-pine)
- [Input file description](#input-file-description)

### Requirements and Setup
The following tools and dependencies are required to run the tool-
1. Install [Cytoscape](https://cytoscape.org/download.html)\
   ![Cytoscape](Image/cytoscape.jpg)
1. Install Cytoscape Apps\
   To install apps within Cytoscape navigate to Apps->App Manager on the tab at the top of the Cytoscape screen. Install the following apps:
   - Genemania
     ![Genemania installation](Image/genemania.jpg)
   - ClueGO (Requires license for usage. Once installation is complete, opening the app prompts license registration)
     ![ClueGO installation](Image/cluego.jpg)
     Once installed apps can be opened by navigating to Apps-> [App Name] on the tab at the top of the Cytoscape screen.
1. Species installation within apps\
   PINE currently supports human, mouse and rat analysis. These species datasets must be installed within the following apps:
   - Genemania
     ![Genemania species installation](Image/genemania-species-install.jpg)
   - ClueGO (by default human and mouse datasets are installed; all other datasets for supported organisms must be installed manually)
     ![ClueGO species installation](Image/cluego-species-install.jpg)
1. Download and Install PINE.exe
   - Download Pine.Setup.zip file from the the latest [release](https://github.com/Niveda-S/PINE/releases) and extract. Click on the .exe file and follow installation instructions

### Using PINE
1. **Setting up PINE:** when you first launch PINE, it will search your PC for the latest Cytoscape executable and ClueGO configuration directory. If it cannot find them, then you will need to manually provide them.  These settings will be saved so they only need to be provided the first time you use PINE.
   ![PINE setup usage](Image/pine-usage-setup-1.png)
1. **Running an analysis:** to begin an analysis, go to the Settings tab.  The required options must be provided before the analysis can be started.
   ![PINE setup usage](Image/pine-usage-settings-1.png)
1. **Analysis options:** the following are the options that can be set to run an analysis.
   - **Input file** - The input file in csv format. See more about the input file in the [input file description section](#input-file-description)
   - **Type** - The type of analysis that will be run. See the [input file description section](#input-file-description) for which types of analysis can be run with which types of input files.
   - **Output directory** - Path to the output directory. A subdirectory for your results will be created within the output directory. The results subdirectory can optionally be named using the input box at the bottom of the settings page. If the subdirectory is not named, then it will be given a name based on the current time.
   - **Species** - Which species are used in your analysis. Currently, only mouse, human and rat are supported. The species must be installed by ClueGO and GeneMANIA in order to use that species. By default ClueGO installs with mouse and human.
   - **Fold change cutoff** - The fold change cutoff for the input. All fold change or log fold change values `>= cutoff` or `<= -cutoff` will be retained.
   - **P-value cutoff** - The p-value cutoff for the input. All p-values or adjusted p-value `<= cutoff` will be retained. If both are provided, adjusted p-value takes precedence.
   - **Outline significant** - Outline statistically significant nodes (p-value or adjusted p-value <= 0.05).
   - **Exclude ambiguity** - Remove all protein, gene and PTM site level ambiguity. (See ambiguity warning in the log).
   - **Database** - Which protein-protein interaction database to use.
   - **Confidence score** - Interaction confidence score for STRING. Highest = 0.9, high = 0.7, medium = 0.4, low = 0.15
   - **# of interactors** - Maximum number of external interactors.
   - **Visualize** - Ontology type. Pathways include REACTOME, KEGG, CLINVAR and Wiki. GO terms for molecular function, biological process and cellular component.
   - **ClueGO grouping** - Network specificity indicating general, representative and specific pathways.
   - **ClueGO p-value** - P-value cutoff for enrichment analysis.
   - **Reference file** - Background reference file containing a list of protein or gene IDs for enrichment analysis in text format.

### Input file description