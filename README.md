# YarraModules-GRASP-Basic

Yarra integrated Matlab implementation of the basic GRASP reconstruction technique. Developed by Li Feng, Ricardo Otazo, and Kai Tobias Block.

This module uses several open-source Matlab packages developed by other contributors, including the mapVBVD library by Philip Ehses, the NUFFT library developed by Jeff Fessler, and the inifile reader developed by Primoz Cermelj. Please give credit to these developers when using the module.


## Reference

This method implements the "GRASP" reconstruction for free-breathing DCE-MRI acquisitions, which has been described in the following publication:

Feng L, Grimm R, Block KT, Chandarana H, Kim S, Xu J, Axel L, Sodickson DK, Otazo R. 
Golden-angle radial sparse parallel MRI: combination of compressed sensing, parallel imaging, and golden-angle radial sampling for fast and flexible dynamic volumetric MRI. 
Magn Reson Med. 2014 Sep;72(3):707-17.

https://www.ncbi.nlm.nih.gov/pubmed/24142845

Please cite this publication when using the module for your work.


## License Information
The Yarra framework is provided free of charge for use in research applications. It must not be used for diagnostic applications. The author takes no responsibility of any kind for the accuracy or integrity of the created data sets. No data created using the framework should be used to make diagnostic decisions.

The Yarra framework, including all of its software components, comes without any warranties. Operation of the software is solely on the user's own risk. The author takes no responsibility for damage or misfunction of any kind caused by the software, including permanent damages to MR systems, temporary unavailability of MR systems, and loss or corruption of research or patient data. This applies to all software components included in the Yarra software package.

The source is released under the GNU General Public License, GPL (http://www.gnu.org/copyleft/gpl.html). The source code is provided without warranties of any kind.

## More Information
More information can be found on the project website http://yarraframework.com

## Notes on Modifications
The GRASP reconstruction code was modified to work with multiple measurement TWIX files.
First the code loads the "mode_file" to overwrite the default settings. Any change on the reconstruction parameters should be done on the "mode_file".
The reconstruction paramters are:
    1) spokesperframe
    2) lambda
    3) slicefrom
    4) sliceto
    5) maxthreads
    6) applykzfilter
All of these parameters can be manualy changed on the mode file "Grasp_basic.mode" found on the "modes" folder.
Parameters can be instroduced under the "[GRASP]" tag.

The TWIX file is read using the "mapVBVD" function found on the "code" folder.
mapVBVD organizes the scan data as:
    1. Columns
    2. Channels/Coils
    3. Lines
    4. Partitions
    5. Slices
    6. Averages
    7. (Cardiac-) Phases
    8. Contrasts/Echoes
    9. Measurements
    10. Sets
    11. Segments
    12. Ida
    13. Idb
    14. Idc
    15. Idd
    16. Ide
For data arrays with more than 4 dimensions all the information corresponding to dimension higher than 4 are grouped together into dimension 5. Where information on dimension 5 is squeezed together.
The TWIX files used in this project contained no information on tags 5-16, with expection of tag 9 corresponding tp the number of measurements. Therefore, the 5 dimension if the data array corresponded to the number of measurements performed.
A for-loop was introduce to iterate the data of the data array using the 5th dimension.
Furthermore, the TA (total scan time) is the combined acquisition time of all measurements; therefore, TA should be adjusted for the corresponding number of scans within the raw data (line 122).
For a more detail description of "mapVBVD" refere to the README file on the code's corresponding folder.

Given that many of the dicom headers are lost due to deidentification and reconstruction of the volunteer data special care is needed when saving the reconstructed scans.
Within the "rave_recon.m" file the variable "output_path" is created, where the output folder should be unique to the desired recontructed scan. Here a folder with name "n_Spokes" will be created by the code for a recontruction using "n" numbers of spokes per frame.

Due to problems with Artemis HPC the parallel processing option was not used, which limits the code's running time. If this issue is resolved uncomment this section (lines 193-197).

Finally, if there is a need to recontruct the data using a list of spokes per frames uncomment the sections on "yarra_GRASP_basic_FN.m" and "rave_recon.m" marked with a "(1)" and comment those marked with a "(2)" as instructed in the corresponding files.

## Artemis HPC (rave.pbs)
To run the code using Artemis HPC the "rave.pbs" should be use as the instruction file.
To enable email notifications to track the job's updates introduce the decired email address under the "PBS -M" tag.
MatLab 2020a was used to run reconstruction code with 32 cpus corres and 123 gb of memory.
For further inquiries on refer to the following link:
https://sydneyuni.atlassian.net/wiki/spaces/RC/pages/185729027/Getting+Started+with+Artemis+HPC