MIST
----

# New in Version 1.5

MiST now automatically scales the specificity score to dataset size (number of independent baits). See `docs/ScalingSpecificity.pdf` for a discussion of this change. This means MiST thresholds and weights should be closer to those used in the original study, but also scores are sensitive to MiST version. *Don't change versions mid-study!*

To disable this specificity score scaling, change your config yml file to include:

```
  standardize_specificity : FALSE
```

See `tests/small/mist_small_test.yml` for example config yml file.

# 1. Necessary Resources

## 1.1 Hardware

Workstation running any current OS, Unix environment recommended

## 1.2 Software

- R package (http://www.r-project.org)
- R packages: getopt, optparse, reshape2, pheatmap, RcolorBrewer, ggplot2, MESS, yaml
- MiST source code (https://github.com/everschueren/MiST)
- Git (optional) (http://git-scm.com)

# 2. Installation

- Download the MiST package as a .zip archive from the public GitHub repository by clicking on the “Download ZIP” button on the bottom right, unzip the files and
- move the directory to a permanent location.
+ Alternatively, you can check out the MiST package through Git as follows:
+ git clone https://github.com/everschueren/MiST.git MiST
- The MiST pipeline is designed to run from a terminal using R. This requires the user to have executable permissions. To set these permissions in a Unix environment, navigate in the terminal to the MiST directory, hereafter referred to as the $INSTALL_DIR, then type: sudo chmod -R 775 *

# 3. Test 

To test your installation, navigate to the `./tests/small` directory within the MIST downloads and run MIST there. On MacOSX and other unix terminals, commands would look like this:

```
cd mist/tests/small
Rscript ../../main.R  -c mist_small_test.yml	
``` 

Output will appear within a subdirectory of your current working directory at `processed/preprocessed_NoC_MAT_MIST.txt`


# 4. References

- Jäger, S., Cimermancic, P., Gulbahce, N. et al. Global landscape of HIV–human protein complexes. Nature 481, 365–370 (2012). https://doi.org/10.1038/nature10719 (See [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fnature10719/MediaObjects/41586_2012_BFnature10719_MOESM288_ESM.pdf), page 46)
- Verschueren, E., Von Dollen, J., Cimermancic, P., Gulbahce, N., Sali, A., & Krogan, N. J. (2015). Scoring Large-Scale Affinity Purification Mass Spectrometry Datasets with MiST. Current protocols in bioinformatics, 49, 8.19.1–8.19.16. https://doi.org/10.1002/0471250953.bi0819s49. (May be [freely available at PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4378866/))
