copy VPRO files
ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("cp $f VPRO/LSQ/vpro_$_\n");' > run_files/copy_vpro.bash

