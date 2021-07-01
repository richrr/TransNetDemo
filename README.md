# TransNetDemo
R package for an introduction to transkingdom network analysis


# July 1 2021:
Since ProNet is no longer available on CRAN, please follow the below steps to install TransNetDemo:

- Download zipped code folder from https://github.com/richrr/TransNetDemo and unzip, place in the user’s local home folder
- Point to the file “DESCRIPTION” and Ctrl + ”open with” a text editor
- Delete the line with ProNet and save the file “DESCRIPTION”
- In R 4.0, install TransNetDemo with devtools:
> devtools::install("/path to TransNetDemo-master/")

You can comment the line "Identify_subnetworks(outNetwork)" and use Cytoscape plugin to identify them using MCODE or InfoMap.


Comment:
- Remotes: AckerDWM/gg3D
- This line in DESCRIPTION tells R that gg3D is on github in AckerDWM/gg3D
- This tells you how to include data in your package: https://r-pkgs.org/data.html


