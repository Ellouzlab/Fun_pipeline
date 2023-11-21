# Fun_pipeline
Automated pipeline for annotating fungi genomes built with funannotate.

Prior to running change fun_exports.sh to match the paths to program and database in your system. Then run:

`source fun_exports.sh`

To rerun from a certain step, delete the file/folder that step produces. For example, to rerun functional annotation, delete Pseudogenus_specicus_strain_XX_predict_folder/annotate_results folder. Halted runs can be continued from checkpoints in this manner.
