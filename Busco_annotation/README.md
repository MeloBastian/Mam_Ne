# **/Busco_annotation** (total size = 10K)

We annotated the 198 genomes using the BUSCO software. 
For computational facilities, we split the BUSCO database (around 9200 genes profiles, orthologous to mammalia) in 27 splits and merge the results in a second step.
      
- ### **/script**

  - **/localbusco.slurm** (input: <genome>, <busco_split>; output: <busco_output> ): Run the Busco analysis.\
    Busco create a lot of temporary files that that habe to be delete upon completion of the Busco run, in order to save file space when parallelising the analysis (implemented by the cleanup function in the script)
       
- ### **/data**
  
  - **/Enard_acc_199sp.csv** : a table with the genome accession id for the 198 species.
    
  - Busco raw output (available on request). Correspond to 197*27 Busco run (=103Go)
