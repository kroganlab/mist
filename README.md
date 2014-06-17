MIST
----

1. Necessary Resources

1.1 Hardware

Workstation running any current OS, Unix environment recommended

1.2 Software

-R package (http://www.r-project.org)
-R packages: getopt, optparse, reshape2, pheatmap, RcolorBrewer, ggplot2, MESS, yaml
-MiST source code (https://github.com/everschueren/MiST)
-Git (optional) (http://git-scm.com)

2. Installation

- Download the MiST package as a .zip archive from the public GitHub repository by clicking on the “Download ZIP” button on the bottom right, unzip the files and
- move the directory to a permanent location.
+ Alternatively, you can check out the MiST package through Git as follows:
+ git clone https://github.com/everschueren/MiST.git MiST
- The MiST pipeline is designed to run from a terminal using R. This requires the user to have executable permissions. To set these permissions in a Unix environment, navigate in the terminal to the MiST directory, hereafter referred to as the $INSTALL_DIR, then type: sudo chmod -R 775 *