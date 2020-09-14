# PINE Installation and Usage

## Table of Contents
- [What is PINE](#what-is-pine)
- [Requirements and Setup](#requirements-and-setup)
- [Example usage](#example-usage)
- [Using PINE GUI](#using-pine-gui)
- [Ontology and Pathway Term Status](#ontology-and-pathway-term-status)
- [Customized Styling in Cytoscape](#customized-styling-in-cytoscape)
- [Handling Ambiguity and PTMs](#handling-ambiguity-and-ptms)
- [Using PINE command line](#using-pine-command-line)
- [Input file description](#input-file-description)
- [Output directory description](#output-directory-description)
- [Cite us](#cite-us)
- [Support](#support)
- [Release notes](#release-notes)
- [License](#license)

## What is PINE
PINE (**P**rotein **I**nteraction **N**etwork **E**xtractor) is a tool for visualizing protein-protein interactions.  PINE is provided in two forms: a command line version built in Python and a graphical user interface (GUI).

## Requirements and Setup
The following tools and dependencies are required to run PINE:
- [Cytoscape](https://cytoscape.org/download.html) (version 3.7 and above)
- [Genemania](http://apps.cytoscape.org/apps/genemania) (version 3.5 and above)
- [ClueGO](http://apps.cytoscape.org/apps/cluego) (version 2.5 and above)

PINE has been tested on Windows 10 and is only available for Windows OS.  We recommend at least 8GB of memory to run PINE with Cytoscape.

### Install [Cytoscape](https://cytoscape.org/download.html)
![Cytoscape](Image/cytoscape.jpg)

### Install Cytoscape Apps
To install apps, open Cytoscape and navigate to Apps->App Manager on the tab at the top of the Cytoscape window. Install the following apps:
#### [Genemania](http://apps.cytoscape.org/apps/genemania)
![Genemania installation](Image/genemania.jpg)
#### [ClueGO](http://apps.cytoscape.org/apps/cluego)
ClueGO requires license for usage. Once installation is complete, opening the app prompts license registration

![ClueGO installation](Image/cluego.jpg)

Once installed apps can be opened by navigating to Apps-> [App Name] on the tab at the top of the Cytoscape screen.

### Species installation within apps
PINE currently supports human, mouse, rat, arabidopsis thaliana, bovine, dog, zebrafish, e.coli, chicken, rabbit, sheep, yeast and pig analysis. These species datasets must be installed within the following apps:
#### Genemania

If you've never used Genemania before, navigate to Apps -> GeneMANIA -> Local Search...

![Genemania species installation](Image/genemania-initial-install-1.png)

Then select the latest dataset

![Genemania species installation](Image/genemania-initial-install-2.png)

Then select which species you would like to install

![Genemania species installation](Image/genemania-initial-install-3.png)

If you've used Genemania before, then open the App within Cytoscape and click on Install Data...

![Genemania species installation](Image/genemania-species-install.jpg)
####  ClueGO
By default human and mouse datasets are installed; other datasets must be installed manually.
If you installed the other datasets in ClueGO, but it still cannot used in PINE, then there may have been a problem during installation of the dataset.
Try reinstalling the dataset to fix the issue.

![ClueGO species installation](Image/cluego-species-install.png)

### Download PINE
Download and run pine-setup.exe from the latest [release](https://github.com/csmc-vaneykjlab/pine/releases/latest/download/pine-setup.exe).

### Windows Defender Smartscreen
You may encounter the following error when running the PINE setup:

![Windows Defender Smartscreen](Image/smartscreen-1.png)

To allow PINE to install, right click on the PINE icon within your Downloads folder and select Properties (if you clicked "Run" instaed of "Download", you will have to [download](https://github.com/csmc-vaneykjlab/pine/releases/latest/download/pine-setup.exe) the file again and save it in a folder on your system):

![Windows Defender Smartscreen](Image/smartscreen-2.png)

At the bottom of the Properties in the Security section, click the checkbox labeled "Unblock" then click "Apply":

![Windows Defender Smartscreen](Image/smartscreen-3.png)

Windows Defender should now allow you to install PINE.

### Install PINE
Open the pine-setup.exe file and follow installation instructions.

Choose if you want to install for a single user account or for all users on the system:

![Installation](Image/install-1.png)

Select the installation directory:

![Installation](Image/install-2.png)

After installation completes you can run PINE:

![Installation](Image/install-3.png)

## Example usage
Please refer to our single fold change [example](examples/single%20fold%20change%20PTM) for a walkthrough on how to use PINE.

## Using PINE GUI

### Setting up PINE
When you first launch PINE, it will search your PC for the latest Cytoscape executable and ClueGO configuration directory. If it cannot find them, then you will need to manually provide them.  These settings will be saved so they only need to be provided the first time you use PINE.

![PINE setup usage](Image/pine-usage-setup-1.png)

### Running an analysis
To begin an analysis, go to the Settings tab.  All required options must be provided before the analysis can be started (see below for an explanation of each option).  Click **Start** to run the analysis.  Click **Load session** to load a previous session.

![PINE setup usage](Image/pine-usage-settings-1.png)

#### Analysis options
The following are the options that can be set to run an analysis.
- **Input file** (required): The input file in csv format. See more about the input file in the [input file description section](#input-file-description)
- **Type** (required): The type of analysis that will be run. See the [input file description section](#input-file-description) for which types of analysis can be run with which types of input files.
- **Output directory** (required): Path to the output directory. A subdirectory for your results will be created within the output directory. The results subdirectory can optionally be named using the input box at the bottom of the settings page. If the subdirectory is not named, then it will be given a name based on the current time.
- **Species** (required): Which species are used in your analysis. Currently, only mouse, human and rat are supported. The species must be installed by ClueGO and GeneMANIA in order to use that species. By default ClueGO installs with mouse and human.
- **Fold change cutoff**: The fold change cutoff for the input. All fold change or log fold change values `>= cutoff` or `<= -cutoff` will be retained.
- **P-value cutoff**: The p-value cutoff for the input. All p-values or adjusted p-value `<= cutoff` will be retained. If both are provided, adjusted p-value takes precedence.
- **Outline significant**: Outline statistically significant nodes (p-value or adjusted p-value <= 0.05).
- **Exclude ambiguity**: Remove all protein, gene and PTM site level ambiguity. (See ambiguity warning in the log).
- **Database**: Which protein-protein interaction database to use.
- **Confidence score**: Interaction confidence score for STRING. Highest = 0.9, high = 0.7, medium = 0.4, low = 0.15
- **# of interactors**: Maximum number of external interactors.
- **Visualize**: Ontology type. Pathways include REACTOME, KEGG, CLINVAR and Wiki. GO terms for molecular function, biological process and cellular component.
- **ClueGO grouping**: Network specificity indicating general, representative and specific pathways.
- **ClueGO p-value**: P-value cutoff for enrichment analysis.
- **Reference file**: Background reference file containing a list of protein or gene IDs for enrichment analysis in text format.

#### PTM analysis options
 If a PTM type analysis is selected, then three more required options will be needed to run the analysis.
- **Fasta file** (required): Fasta file for finding locations of PTMs within the protein.
- **Enzyme** (required): Digestion enzyme used.
- **Modifications** (required): Comma separated list of modifications of interest.

### Log
While PINE is running, output is written to the log. A copy of the log will also be saved to PINE.log within the output subdirectory. The analysis can be cancelled at anytime from this tab.

![PINE log example](Image/pine-usage-log-1.png)  
When the PINE analysis is complete, you can view the interaction network by going to the Network tab and selecting the Interaction Network.

![Interaction network example](Image/initial-network.png)

### Pathway selection
After a PINE analysis successfully completes, the pathway selection tab will load which shows all the pathways and GO terms found from ClueGO analysis. Pathways and GO terms can be selected for reanalysis on the subset of genes found within these terms. The reanalysis can be given a custom name. If a name is not given, then it will automatically named based on the current time. After selecting one or more terms or pathways, click "Reanalyze" to begin the reanalysis.

From this tab, you can also open the results folder ([see here for results description](#output-directory-description)) and the Cytoscape file of the most recent analysis.

![Pathway analysis](Image/pine-usage-pathway-selection-1.png)

The columns in the pathway selection page are:
- **GO term**: Annotation or pathway term for a query
- **P-value**: Significance of the GO term
- **Adj.p-value**: Corrected significance of the GO term
- **% genes**: Percent of genes in a ClueGO cluster associated with the GO term
- **# of genes**: Number of genes associated with the GO term
- **Status**: Depiction of up-regulation, down-regulation or no change status of a GO term determined on the basis of status of a majority of associated genes. If a majority of genes have fold changes > 1, then the GO term is said to be overall up-regulated. If a majority of enes have fold changes < 1, then the GO term is said to be overall down-regulated. If there is no majority of genes with either fold change > 1 or < 1, then the GO term is said to have no change.

After reanalysis is complete, there will be a new interaction network which contains only the genes from the selected pathways and terms and an ontology network shows which genes are included in the selected pathways and terms.

![Pathway analysis ontology network](Image/pathway-selection-network.png)

Colors on the annotation node indicate the following:
- **Dark blue annotation nodes**: Nodes connected to 60% or more up-regulated sites
- **Orange annotation nodes**: Nodes connected to 60% or more down-regulated sites
- **Grey annotation nodes**: Nodes connected to approximately the same number of up-regulated and down-regulated sites  

**NOTE**: Annotation nodes are fixed colors and do not indicate the degree of connectivity to up-regulated or down-regulated sites  

### Legend  
Shown below is the legend for interpreting visualizations through PINE:
![Legend](Image/legend.png)

## Ontology and Pathway Term Status
In the case of single-fc type analysis, a column labelled "Status" is available in the Pathway Selection page that depicts up-regulation, down-regulation or no change status of a GO term determined on the basis of status of a majority of associated genes. If a majority of genes have fold changes > 1, then the GO term is said to be overall up-regulated. If a majority of enes have fold changes < 1, then the GO term is said to be overall down-regulated. If there is no majority of genes with either fold change > 1 or < 1, then the GO term is said to have no change.

![Pathway analysis singlefc](Image/pine-usage-pathway-selection-1.png)

In the case of multiple-fc type analysis, a drop-down of all comparison groups appears on the top of the Pathway Selection page. Based on the label selected, a column with the comparison name is updated that depicts up-regulation, down-regulation or no change status of a GO term determined on the basis of status of a majority of associated genes. Switching comparison groups in the dropdown updates the column to the selected comparison name.

![Pathway analysis multifc](Image/multifc-pathway-selection.png)

## Customized Styling in Cytoscape

### Layout
It is recommended to install the Cytoscape App yFiles Layout Algorithms to improve your layouts after PINE has finished building your network. To install navigate to Apps->App Manager on the tab at the top of the Cytoscape screen.
Appropriate yFiles layout can be chosen by selecting from options available in the Layout tab at the top of the Cytoscape screen. For example, figures 4 and 5 in the manuscript have been constructed using yFiles Organic Layout. Figures 6 and 7 have been constructed using yFiles Heirarchic Layout. 

### Font
To change font, navigate to the Style Tab on the Control Panel on the left side of the Cytoscape screen. 'Label Font Size' is the option to set node font size. Font has been set by PINE based on length of node text. To change, the range of minumum and maximum font must be modified.
Open up the node font size section by clicking on the dropdown arrow. Double-click on the current mapping graph. This opens up a new tab called Continuous Mapping Editor. Move the node label font size pointer up or down to increase or decrease font size respectively.
![Font Size](Image/FontSize.png)

Alternatively, to change font sizes of specific node labels, click to select the required node and select the bypass box for Label Font Size on the Style tab. This opens up a new tab where node label font size for the selected node can be set.
![Bypass](Image/BypassFont.png)

### Colors
**Color Gradient for single fold change type networks**  
To modify color of fold change greadient for the nodes, navigate to the Style Tab on the Control Panel on the left side of the Cytoscape screen. 'Fill Color' controls node colors.   
Open up the fill color section by clicking on the dropdown arrow. Double-click on the current mapping gradient. This opens up a new tab called Continuous Mapping Editor. Set a new gradient by double-clicking on the arrow buttons for fold changes that are up-regulated (>0) and down regulated(<0) and setting new colors from a color palette.  
![SinglefcColor](Image/SinglefcColor.png)

**Bar and Pie Charts for multiple fold change or category type networks**  
The colors of the bar and pie charts can be modified by navigating to the 'Image/Chart' option on the Style Tab. Click on the currently configured bar chart to open up a Graphics tab. Its options section shows the currently configured colors. Double-clicking on the colors opens up a color palette allowing you to make any changes.  
![BarColor](Image/BarChartColor.png)  

Additionally, for bar charts, the option is provided to include:  
- **Value Labels** in the bar showing the numeric value that represents height of the bar. Select 'Show Value Labels' and choose the  column --none-- to enable this option  
- **Domain Labels** showing bar labels. Select 'Show Domain Axis' and choose the PINE generated column pine_domain_label to enable this option  
Other options include ability to change font size of labels and values, showing axis lines and altering width of axis.
![BarOptions](Image/BarChartOptions.png)

This labels bar charts appropriately as below:  
![Ex](Image/ExBarChart.png)  

NOTE: In order to change the numeric values depicted in the bar charts, make the required changes to the Cytoscape Fold Change columns of the appropriate label.

Further, the position of the horizontal axis of the bar charts can also be altered. By default, the lowest and highest points of the bar chart are determined by the minimum and maximum values of the input dataset. This can be modified to adjust the bar charts to the centre of the node in cases where minimum or maximum values are 0 indicating that the bar charts either start or end at the bottom of node respectively. Decreasing the minimum shifts the horizontal axis to the top and increasing the maximum moves the axis to the bottom:
![min-and-max](Image/min-and-max-bar-chart.jpg)
Uncheck and recheck the 'Automatic Range' option to reset minimum and maximum values

Note: 
1. Changing the minimum and maximum values too much from the minimum and maximum values of the dataset could impact the height of the bars and hamper your ability to differentite between values represented by the bar.
2. Care must be taken in case of repeated checking and unchecking of the 'Automatic Range' option as the minimum and maximum values automatically gets set to -100 and +100 respectively.
3. Addition/Deletion of columns from the bar chart comparison is possible by moving the columns of interest to/from 'Selected Columns' tab. In this case, minimum and maximum values may automatically gets set to -100 and +100 respectively and must be reset appropriately to view the bar chart.

## Handling Ambiguity and PTMs
### Ambiguity Resolution
| Case |	Category |	Step |	Level |	Description |	Resolution |
| ------ | ------ | ------ | ------ | ---------------- | ---------------- |
| Ambiguous sites	| Ambiguity	| Preprocessing |	PTM Site |	Multiple peptides that represent same PTM type on the same site of a protein |	One representative peptide picked |
| Peptide Mapping |	Ambiguity |	Preprocessing |	PTM Site |	Peptide maps to multiple regions in FASTA |	First mapping site picked |
| Isoforms |	Ambiguity |	Uniprot Mapping |	Protein |	Isoforms in query map to same primary gene |	One represenatative isoform picked, canonical protein, if it exists gets preference |
| Primary Gene Mapping |	Ambiguity |	Uniprot Mapping |	Protein	| Query ProteinID maps to multiple primary genes in Uniprot |	First mapping primary gene picked |
| Duplicate Primary Gene	| Ambiguity	| Uniprot Mapping	| Protein	| Query ProteinIDs map to same primary gene	| One representative protein ID picked, reviewed protein, if it exists gets preference |
| Duplicate Protein Mapping | Ambiguity | Uniprot Mapping | Protein | Multiple query ProteinIDs mapping to single Uniprot ID | Drop all obsolete IDs, retain active IDs |
| Peptide Unmapped |	Discard	| Preprocessing	| Peptide |	Query peptides do not map in FASTA	| Drop peptides with no peptide mapping |
| Duplicate Peptide | Discard | Preprocessing | Peptide | Duplicate query peptides across multiple UniprotIDs | Drop duplicate peptides |
| Site Not Available |	Discard	| Preprocessing |	Peptide |	No modifications of interest found in peptide	| Drop all peptides having no modification site |
| Duplicate query	| Discard |	Preprocessing |	Protein |	Duplicate fields in input |	Drop duplicates |
| Invalid FC/Pval	| Discard	| Preprocessing	| Protein |	Query contains non-numeric fold change and p-values |	Drop all queries with non-numeric FC/Pval | 
| FC/Pval cutoff |	Discard |	Preprocessing |	Protein |	Query does not meet fold change and p-value cutoffs |	Drop all queries with cutoff not met |
| ProteinID unmapped |	Discard	| Uniprot Mapping |	Protein	| Query ProteinID not mapped in Uniprot	| Drop queries with no Uniprot mapping |
| Protein Mapping |	Discard |	Uniprot Mapping |	Protein |	Query ProteinID maps to multiple UniprotIDs |	Drop all queries because ID is obsolete |
| Primary Gene Not Available |	Discard |	Uniprot Mapping |	Protein	| Query ProteinID does not have a primary gene in Uniprot	| Drop queries with no primary genes |
| Gene Unmapped in Interaction Databases	| Discard	| Interaction Retrieval	| Gene	| Genes not mapped in String and GeneMANIA	| Drop all unmapped genes |
| Interactions Not available	| Discard	| Interaction Retrieval	| Gene	| Genes having no interactions in String and GeneMANIA	| Drop all genes with no interactions |
| Invalid Interaction Category	| Discard	| Interaction Retrieval	| Gene	| Query genes not categorized as primary interactors in STRING and GeneMANIA	| Drop all query genes with category other than primary interactor |
| Gene Unmapped in ClueGO	| Discard	| Enrichment	| Gene	| Genes not mapped in ClueGO	| Drop all unmapped genes |
| Invalid Mapping Category	| Discard	| Enrichment	| Gene	| Genes not categorized as primary in ClueGO mapping	| Drop all query genes with category other than primary |

### PTM Handling  
For ambiguous sites, scoring is performed for representative site selection by excluding sites with oxidized methionines, missed cleavages, and ragged ends.  
- Missed Cleavages  
![Missed-Cleavage](Image/Missed-Cleavage.JPG)  
Number of missed cleavages are calculated for all peptides of a single site, with the peptide having the least number of missed cleavages showing higher score.  

- Other Modifications  
![Other-Mods](Image/Other-Mods.JPG)  
Number of modifications other than the modification of interest are calculated for all peptides of a single site, with the peptide having the least number of other modifications showing higher score.  

Note: Alkylated Cys are not considered as other modifications. Sequences containing alkylated Cys are always preferred over sequences containing non-alkylated Cys during ambiguity resolution.

### PTM Naming Convention
PTM sites are represented in the interaction and ontology network by Amino acid modified followed by PTM type in curly brackets followed by PTM site as shown above (E.g. S{+80}25). PTM type is denoted based on PTM identifier present within brackets in the input data (e.g. Modification mass [+80] or Unimod accession (Unimod:21) or free text {Phos} etc.). This information proves to be useful in case of occurrence of multiple PTMs on a single amino acid in order to differentiate between PTMs using the PTM type. But, in other cases, this can be optionally turned off by the user by switching node label to column 'substitute name' within Cytoscape->Control Panel->Style->Label (e.g. S{+80}25 will be denoted as S25) as shown below:  
![substitute-name](Image/substitute-name.png)  

NOTE: Ensure that the 'Mapping Type' of the node label remains Passthrough Mapping:
![passthrough](Image/passthrough.png)

## Using PINE command line
### Requirements
- All the requirements listed [above](#requirements-and-setup)
- Python3

### Installation
`pip install git+https://github.com/csmc-vaneykjlab/pine.git#egg=pine`

It is strongly recommended you use a virtual environment for installing PINE.

After installing, PINE can be run using `python3 -m pine.pine --help`

#### Windows issues
Installing the python-igraph dependency may have issues on Windows.
If you encounter an error related to python-igraph, you can install a binary of the package from here: https://www.lfd.uci.edu/~gohlke/pythonlibs/#python-igraph.
After downloading the package, run `pip install path\to\wheel\python_igraph‑0.7.1.post6‑cp38‑cp38‑win_amd64.whl` then rerun the install command from above.

### Usage
```
python3 -m pine.pine -i input.csv -o output_dir -c cluego_out.txt -t input_type -s species -m cluego_map_file.gz --cytoscape-executable path_to_exe
```

### Command line parameters
| Parameter | Description |
| --------- | ----------- |
| -i, --in | input file in csv format with the following headers as applicable: ProteinID, GeneID, FC, pvalue, adj.pvalue, Label, Category, Peptide |
| -o, --output | path to output directory |
| -t, --type | analysis type [Allowed: noFC, singleFC, multiFC, category, singlefc-ptm, multifc-ptm] |
| -s, --species | species [Allowed: human, mouse, rat, arabidopsis, bovine, dog, zebrafish, e. coli, chicken, rabbit, sheep, yeast, pig] |
| -x, --enzyme | (required if singlefc-ptm or multifc-ptm) enzyme name [Allowed: Trypsin, Trypsin_p, Lys_n, Asp_n, Arg_c, Chymotrypsin, Lys_c] |
| -d, --mods | (required if singlefc-ptm or multifc-ptm) comma separated list of modifications of interest [Example: S,T,Y or K(Unimod:1) or S[+80]] |
| -b, --fastafile | (required if singlefc-ptm or multifc-ptm) path to fasta file |
| -m, --mapping | path to cluego mapping file compressed in .gz format |
| -e, --cytoscape-executable | the path to the Cytoscape executable |
| -f, --fccutoff | (optional) fold change cutoff for input [Default: abs(FC) >= 0.0] |
| -p, --pvalcutoff | (optional) pvalue cutoff for input [Default: pval > 1.0] |
| -n, --significant | (optional) outline statistically significant nodes, i.e pval>0.0 |
| -k, --exclude-ambiguity | (optional) exclude ambigious genes and sites |
| -u, --run | (optional) interaction databases [Allowed: string, genemania, both; Default: both] |
| -r, --score | (optional) interaction confidence score for string [Default:0.4, Range 0-1] |
| -l, --limit | (optional) maximum number of external interactors [Default:0, Range:0-100] |
| -z, --visualize | (optional) ontology type [Allowed: biological process, cellular component, molecular function, pathways, all; Default: pathways].  Pathways include REACTOME, KEGG, CLINVAR, CORUM and Wiki. |
| -g, --grouping | (optional) network specificity indicating general, representative and specific pathways [Allowed: global, medium, detailed; Default: medium] |
| -y, --cluegopval | (optional) pvalue cutoff for enrichment analysis [Default: 0.05] |
| -h, --referencepath | (optional) path to background reference file for enrichment |
| -a, --inputcluego | (optional) filtered cluego file with ontology terms of interest |

## Input file description
All input files must be in CSV (comma separated value) format.  All column names are case-insensitive.

**No fold change**

| Column | Input column name |
| ------ | ---------------- |
| Uniprot ID or Gene ID | `proteinid` or `gene name` |

**Single fold change**

| Column | Input column name |
| ------ | ---------------- |
| Uniprot ID or Gene ID | `proteinid` or `gene name` |
| Fold change | `fc` |
| P.value (opt) | `pvalue` or `adj.pvalue` or `fdr` | 

**Multi fold change**

| Column | Input column name |
| ------ | ---------------- |
| Uniprot ID or Gene ID | `proteinid` or `gene name` |
| Fold change | `fc` |
| Label | `label` |
| P.value (opt) | `pvalue` or `adj.pvalue` or `fdr` | 

**Category**

| Column | Input column name |
| ------ | ---------------- |
| Uniprot ID or Gene ID | `proteinid` or `gene name` |
| Category | `category` |

**Single fold change PTM**

| Column | Input column name |
| ------ | ---------------- |
| Uniprot ID | `proteinid` |
| Peptide sequence | `peptide` |
| Fold change | `fc` |
| P.value (opt) | `pvalue` or `adj.pvalue` or `fdr` | 

**Multi fold change PTM**

| Column | Input column name |
| ------ | ---------------- |
| Uniprot ID | `proteinid` |
| Peptide sequence | `peptide` |
| Fold change | `fc` |
| Label | `label` |
| P.value (opt) | `pvalue` or `adj.pvalue` or `fdr` | 

NOTE: Alkylated cys should have the following naming convention - C[+57], C(Unimod:4) or C[cys] 

## Output directory description

A directory is created in the specified output directory after the analysis completes.  This directory will contain six files:
- **Interactions.csv** - The results for each Protein ID in the analysis.
- **PINE.cluego.txt** - Pathways and GO terms found to be significant by ClueGO.
- **PINE.cys** - The Cytoscape file containing the interaction network.
- **PINE.log** - Log of analysis.
- **settings-gui.json** - Settings from the GUI that were provided to the PINE Python script.  This file should not be modified because it is used to retrieve the settings when the session is reloaded.
- **timestamp.json** - The time when the analysis was run.  This file should not be modified.

## Cite us
Niveda Sundararaman, James Go, Aaron E. Robinson, José M. Mato, Shelly C. Lu, Jennifer E. Van Eyk, Vidya Venkatraman. "PINE: An Automation Tool to Extract & Visualize Protein-Centric Functional Networks". Journal of the American Society for Mass Spectrometry. (2020) https://pubs.acs.org/doi/pdf/10.1021/jasms.0c00032

## Support
If you have any questions about PINE, please contact us at GroupHeartBioinformaticsSupport@cshs.org.

## Release notes
### Version 2.0.2
- Fixed directory selection for new Electron version.
- Show version number on setup tab in the GUI.

## License
See the [LICENSE](https://github.com/csmc-vaneykjlab/pine/blob/master/LICENSE) file for license rights and limitations (Apache 2.0).
