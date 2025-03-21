- /Busco_annotation 10K
We annotate the genomes using Busco. For computational facilities, we split Busco database (around 9000 genes profiles) in 27 splits and merge the results later.
        -- /script
                -- /localbusco.slurm (input: <genome>, <busco_split>; output: <busco_output> ): Run the Busco analysis
        -? /data
                -- /Enard_acc_199sp.csv : a table with genome accession number
                -? Busco raw output (?). Correspond to 197*27 Busco run (=103Go)
