#MIMET (MIcrobial METabolites) v0.0.1

Background
------

Creating metabolic models from whole genome sequencing projects is routine these days. However, to be able to 
predict microbial metabolites from 16S rRNA sequencing data is something novel. With this intention in mind, we have
developed a pipeline that relying upon the output of [PICRUSt software] (http://picrust.github.io/picrust/) is able to
predict microbial metabolites from 16S data.


Required Packages
------

**Python:**

- [PICRUSt 1.0.0](http://picrust.github.io/picrust/install.html#install)
- [NumPy 1.8.0](http://www.scipy.org/scipylib/download.html)

How to use
------

To see help with pipeline, simply type the following in the terminal:

```python microbial_metabolites.py -h```

This script requires all other scripts to be present in the same working directory. The script also assumes that the input file is also in the working directory. 
