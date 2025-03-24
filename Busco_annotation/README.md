# **/Busco_annotation** (total size = 10K)

We annotate 198 genomes using the BUSCO software. 
For computational facilities, we split the BUSCO database (around 9200 genes profiles, orthologous to mammalia) in 27 splits and merge the results later.
      
- ### **/script**

  - **/localbusco.slurm** (input: <genome>, <busco_split>; output: <busco_output> ): Run the Busco analysis.\
    Busco create a lot of temporary files that you have to delete to preserve your computer memory when you parallelise the analysis (cleanup function in the script)
       
- ### **/data** (_should we provide the busco output ? very long to run but also very heavy to storage_)
  
  - **/Enard_acc_199sp.csv** : a table with the genome accession id for the 198 species.
    
  - Busco raw output (?). Correspond to 197*27 Busco run (=103Go)
