# Single fold change PTM example
This guide will walk you through an example PINE analysis using for single fold changes PTMs.

## Download the example files
You'll need the [single-fold-change-data.csv](https://raw.githubusercontent.com/csmc-vaneykjlab/pine/master/examples/single%20fold%20change%20PTM/single-fold-change-data.csv) and [mouse.fasta](https://raw.githubusercontent.com/csmc-vaneykjlab/pine/master/examples/single%20fold%20change%20PTM/mouse.fasta) to run this example.

## Open PINE and set the run parameters
Make sure PINE was able to find your Cytoscape.exe and ClueGO configuration folder.  You'll need to provide the locations if it can't find them.
- Select the **single-fold-change-data.csv** as your input file.
- Set type to **Single fold change PTM**.
- Select an output directory to write your results.
- Set the species to **mouse**.
- Select the **mouse.fasta** as your fasta file.
- Set the enzyme to **Typsin**.
- Set the modifications to **S,T,Y**.
- Set fold change cutoff to **0**.
- Set p-value cutoff to **1**.
- **Uncheck** outline significant.
- **Uncheck** exclude ambiguity.
- Set database to **both**.
- Set confidence score to **0**.
- Set # of interactors to **0**.
- Set visualize to **all**.
- Set ClueGO grouping to **medium**.
- Set ClueGO p-value to **1**.
- Do not add a reference file.
- Name your analysis directory **single-fold-change-example**.

Your settings should look like the following when you're done.  Click **Start** to run.
![settings](images/settings.png)

## Analysis results
After PINE finishes running, you can go to the **Network** tab in the Cytoscape control panel and select the **Interaction Network**.  This is the single fold change interaction network.  Each yellow square is a gene and the circles represent fold changes for PTMS on the genes' proteins.  Your network should look similar to the image below.  We recommend installing the [yFiles Layout Algorithms plugin](http://apps.cytoscape.org/apps/yfileslayoutalgorithms) and using the yFiles Organic Layout to layout your networks.
![results 1](images/results-1.png)

## Pathway selection
You can now refine your results by selecting different pathways or GO terms to reanalyze with.  Go back to PINE and select a few terms, then click **Reanalyze**.
![pathway selection](images/pathway-selection.png)

## Reanalysis results
After the PINE reanalysis completes, you will have two networks in your Cytoscape session.  The first is a new interaction network containing only genes from the terms you selected.  The second is an ontology network for associations between genes and the selected ontologies.

**New interaction network:**
![results 2](images/results-2.png)

**Ontology network:**
![results 3](images/results-3.png)