#MIMET (MIcrobial METabolites) v0.0.1

Background
------

Creating metabolic models from whole genome sequencing projects is routine these days. However, to be able to 
predict microbial metabolites from 16S rRNA sequencing data alone is something novel. With this intention in mind, we have developed a pipeline which, based upon the output of [PICRUSt software] (http://picrust.github.io/picrust/), is able to predict microbial metabolites from 16S data.


Required Packages
------

**Python:**

- [PICRUSt 1.0.0](http://picrust.github.io/picrust/install.html#install)
- [NumPy 1.8.0](http://www.scipy.org/scipylib/download.html)

How to use
------

To see help with pipeline, simply type the following in the terminal:

```python microbial_metabolites.py -h```

The input file for this pipeline is the output file generated by the [predict_metagenomes.py](http://picrust.github.io/picrust/scripts/predict_metagenomes.html) script of PICRUSt in tab-delimited format. An example file is provided in the src/ folder called ```input_ko_table.tab```.

The pipeline also requires the "reaction" and "compound" KEGG database flat-files to do the relevant mappings and calculations. 
