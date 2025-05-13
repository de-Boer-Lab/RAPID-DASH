The directory contains all the scripts from the publication.

1. get_plasmid_inserts.txt: Script to extract gRNA array sequences from the whole plasmid sequencing reads. Input required is the .fastq file from Nanopore sequencing

2. (clonal_)nanopore Sequencing analysis.ipynb: Jupyter notebook to analyse the Nanopore reads --- plotting gRNA unit distribution, mapping gRNA spacer sequences and calculating the efficiency of RAPID-DASH. Input required is the .pkl file outputted from get_plasmid_insert command.

3. GFPreporterAssay.R: R script to analyse and plot GFP reporter  flow cytometry data. Input required is the .fcs files from flow cytometer.

4. Array-library.R: R script to plot the distribuiton of gRNA spacer sequences within the array library generated from pooled PCA. Input required is the .csv file of gRNA mapping from nanopore sequencing analysis.ipynb.

5. Assembly_efficiency.R: scipt to calculate and plot the efficiency of array assembly from both clonal and bulk plasmid sequencing data. Input required is the number of arrays with different number of gRNA units assembled from nanopore sequencing analysis.ipynb.

### Environment setup for Jupyter notebook

Create a new conda environment using nanopore_analysis.yaml:
```
conda env create --file=nanopore_analysis.yml
conda activate nanopore_analysis
```
Add the environment kernel:
```
python -m ipykernel install --user --name=nanopore_analysis
````
