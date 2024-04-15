  279  srun ./results_postprocess.py --dpath ../../results/vpro_ponly/out/ --outdpath results --run_id vpro_ponly 
  280  ls
  281  ls results
  282  git st
  283  git diff
  284  git st
  285  squeue 
  286  ls ../../results/lsq_ponly/out/
  287  ls ../../results/lsq_ponly/out/|wc
  288  tail -n 100 ../../results/lsq_ponly/log/* | less
  289  rm VPRO/FINAL/*
  290  cp ../../results/lsq_ponly/out/* VPRO/FINAL/
  291  ls VPRO/FINAL/
  292  ls VPRO/FINAL/ |wc
  293  sacct
  294  sacct | less
  295  rm -rf ../../results/vpro_ponly/out/
  296  rm -rf ../../results/vpro_ponly/log/*
  297  ls
  298  ls run_files/
  299  vim run_files/run_vpros_ponly_sbatch.bash 
  300  source  run_files/run_vpros_ponly_sbatch.bash
  301  squeue 
  302  squeue |less
  303  tail -n 100 ../../results/vpro_ponly/log/* | less
  304  ls ../../results/vpro_ponly/out/
  305  sacct | less
  306  srun ./results_postprocess.py --help
  307  srun ./results_postprocess.py --dpath ../../results/vpro_ponly/out/ --outdpath results --run_id vpro_ponly
  308  ls -ltrh results
  309  squeue 
  310  sacct | less
  311  ls ../../results/vpro_ponly/out/
  312  srun ./results_postprocess.py --dpath ../../results/vpro_ponly/out/ --outdpath results --run_id vpro_ponly
  313  ls
  314  ls -ltrh results
  315  git st | less
  316  less VPRO/FINAL/lsq_ponly_ROS_21.json 
  317  ls ../../results/lsq_ponly/log/out_273*| less
  318  grep ROS_21 ../../results/lsq_ponly/log/out_273*
  319  grep ROS ../../results/lsq_ponly/log/out_273*
  320  less ../../results/lsq_ponly/log/out_2732167_21_2732178.out 
  321  grep cost ../../results/lsq_ponly/log/out_2732167_*
  322  grep cost: ../../results/lsq_ponly/log/out_2732167_*
  323  ls ../../results/vpro_ponly/log/out_273*
  324  grep ROS_21 ../../results/vpro_ponly/log/out_273*
  325  tail  ../../results/vpro_ponly/log/out_273* | less
  326  git st
  327  git co -- VPRO/
  328  git st | less
  329  git checkout -- VPRO/FINAL/
  330  git st | less
  331  rm -rf VPRO/FINAL/
  332  git checkout -- VPRO
  333  git st | less
  334  git st 
  335  git diff
  336  exit
  337  squeue 
  338  exit
  339  scancel 
  340  scancel -u $USER
  341  squeue 
  342  squeue -u $USER
  343  exit
  344  squeue 
  345  exit
  346  srun ls
  347  cd RECYCLE_MODEL/
  348  git st
  349  cd recycle_model/
  350  git st
  351  git pull
  352  git st
  353  git pull
  354  cd STORE_model/
  355  ls
  356  cat run_files/run_lsq_ponly_sbatch.bash 
  357  RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
  358  mkdir slurm
  359  cd slurm/
  360   m=MIN
  361  vim $RDIR/run_files/run_lsq_ponly_parallel.bash
  362  cat run_files/run_lsq_ponly_sbatch.bash 
  363  cat $RDIR/run_files/run_lsq_ponly_sbatch.bash 
  364  ls VPRO/X0/*${m}* |  xargs sbatch --partition=hive1d,hive7d,hiveunlim --wrap "echo {}"
  365  vim /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  366  cat $RDIR/run_files/run_lsq_ponly_sbatch.bash 
  367  ls $RDIR/VPRO/X0/*${m}* |  xargs sbatch --partition=hive1d,hive7d,hiveunlim --wrap "echo $m {}"
  368  ls $RDIR/VPRO/X0/*${m}* |  xargs sbatch --partition=hive1d,hive7d,hiveunlim echo $m {}
  369  vim /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  370  ls $RDIR/VPRO/X0/*${m}* |  xargs sbatch --partition=hive1d,hive7d,hiveunlim  $RDIR/run_files/run_lsq_ponly_parallel.bash $m {}
  371  squeue 
  372  ls
  373  less slurm-2784831.out 
  374  ls $RDIR/VPRO/X0/*${m}* |  xargs -n 1 sbatch --partition=hive1d,hive7d,hiveunlim  $RDIR/run_files/run_lsq_ponly_parallel.bash $m {}
  375  ls 
  376  squeue 
  377  ls 
  378  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/VPRO/X0/shgo3_sobol_370_OVERFLOW.json |  xargs -n 1 sbatch --partition=hive1d,hive7d,hiveunlim  $RDIR/run_files/run_lsq_ponly_parallel.bash $m {}
  379  squeue 
  380  sacct
  381  for j in (ls $RDIR/VPRO/X0/*${m}*); do; echo $j ;sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_parallel.bash $m $j"; done
  382  for j in $(ls $RDIR/VPRO/X0/*${m}*); do; echo $j ;sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_parallel.bash $m $j"; done
  383  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; echo sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_parallel.bash $m $j"; done
  384  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_parallel.bash $m $j"; done
  385  sacct
  386  squeue 
  387  ls
  388  less slurm-278*
  389  cat slurm-278* | less
  390  vim /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  391  ls
  392  rm -rf slurm-278*
  393  cat slurm-278* | less
  394  head /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  395  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  396  m=MIN
  397  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_parallel.bash $m $j"; done
  398  squeue 
  399  ls
  400  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  401  ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
  402  baseid=lsq_ponly_${m}
  403  ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
  404  ls $ODIR
  405  mkdir -p $ODIR/log
  406  ls $ODIR
  407  tail slurm-2785* | less
  408  vim  /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  409  rm -rf slurm-2785
  410  rm -rf slurm-2785*
  411  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  412  tmux attach
  413  RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
  414   m=MIN;
  415  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done
  416  tail slurm-* | less
  417  squeue 
  418  squeue  | less
  419  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  420  m=MIXOTROPH
  421  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done
  422  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  423  m=OVERFLOW
  424  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done
  425  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  426  m=ROS
  427  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done
  428  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  429  m=EXOENZYME
  430  for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done
  431  squeue | less
  432  squeue | grep oweiss
  433  squeue | grep oweiss|wc
  434  tail slurm* | less
  435  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  436  tail $ODIR/log/* | less
  437  git st
  438  squeue 
  439  git pull
  440  rm -rf /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly_MIN/log/*
  441  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  442  for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done; ;
  443  for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do for j in $(ls $RDIR/VPRO/X0/*${m}*); do echo $j; sbatch --partition=hive1d,hive7d,hiveunlim --wrap  "$RDIR/run_files/run_lsq_ponly_sbatch.bash $m $j"; done; done
  444  vim  /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  445  squeue 
  446  git diff
  447  git st
  448  ls
  449  srun tail $ODIR/log/* | less
  450  tail $ODIR/log/* | less
  451  cat * | less
  452  squeue 
  453  squeue |less
  454  squeue -u $USER |less
  455  tail $ODIR/log/* | less
  456  tail slurm-* | less
  457  ls $ODIR/log/ | less
  458  ls $ODIR/out | less
  459  sacct| less
  460  ls $ODIR
  461  ls $ODIR/*
  462  /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/
  463  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/
  464  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out
  465  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | less
  466  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly_EXOENZYME/out/* | less
  467  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly_EXOENZYME/log| less
  468  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly_EXOENZYME/log/*| less
  469  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  470  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | less
  471  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly
  472  rm -rf  /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly
  473  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | less
  474  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/log/* | less
  475  tail /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/log/* | less
  476  sacct| less
  477  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | wc
  478  squeue | less
  479  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | wc
  480  squeue -u $USER | less
  481  squeue -u $USER | wc
  482  squeue -u $USER | less
  483  squeue -u $USER | wc
  484  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | wc
  485  squeue -u $USER | wc
  486  squeue -u $USER 
  487  ls -ltr /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | less
  488  tail /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | less
  489  wc -l /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | less
  490  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | head
  491  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print;'
  492  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("mv $f $RDIR/VPRO/LSQ/vpro$_");'
  493  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("mv $f $RDIR/VPRO/LSQ/vpro$_\n");'
  494  cd RECYCLE_MODEL/recycle_model/
  495  git st
  496  git pull
  497  cat /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model/run_files/run_lsq_ponly_sbatch.bash 
  498  RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
  499  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("mv $f $RDIR/VPRO/LSQ/vpro$_\n");'
  500  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("mv $f $env{RDIR}/VPRO/LSQ/vpro$_\n");'
  501  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("mv $f VPRO/LSQ/vpro$_\n");'
  502  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("cp $f VPRO/LSQ/vpro_$_\n");'
  503  mkdir VPRO/LSQ
  504  cd STORE_model/
  505  mkdir VPRO/LSQ
  506  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("cp $f VPRO/LSQ/vpro_$_\n");'
  507  ls /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/lsq_ponly*/out/* | perl -ne 'chomp;$f= $_; s/.*shgo//;s/_sobol_//; print("cp $f VPRO/LSQ/vpro_$_\n");' > run_files/copy_vpro.bash
  508  source run_files/copy_vpro.bash
  509  ls VPRO/LSQ/
  510  ls
  511  git st
  512  ls run_files/hive/
  513  vim run_files/hive/run_vpros_ponly_sbatch.bash 
  514  source  run_files/hive/run_vpros_ponly_sbatch.bash
  515  ls ../../results/vpro_ponly/
  516  ls ../../results/vpro_ponly/out/
  517  rm -rf  ../../results/vpro_ponly/out/vpro_ponly_lsq_ponly_*
  518  ls ../../results/vpro_ponly/out/
  519  ls ../../results/vpro_ponly/log/
  520  squeue -u $USER  
  521  squeue -u $USER  |wc
  522  squeue -u $USER  |less
  523  ls ../../results/vpro_ponly/out/
  524  ls ../../results/vpro_ponly/out/|wc
  525  squeue -u $USER  |less
  526  ls ../../results/vpro_ponly/out/|wc
  527  squeue -u $USER  |wc
  528  ls ../../results/vpro_ponly/out/|wc
  529  squeue -u $USER  |wc
  530  ls ../../results/vpro_ponly/out/|wc
  531  squeue -u $USER  |wc
  532  squeue -u $USER  
  533  ls ../../results/vpro_ponly/out/|wc
  534  vim run_files/README.txt
  535  ls results
  536  rm -rf  results/vpro_ponly_*
  537  srun ./results_postprocess.py --dpath ../../results/vpro_ponly/out/ --out results --run_id vpro_ponly
  538  vim run_files/run_shgo.bash 
  539  ls run_files
  540  vim run_files/run_shgo_ponly_hive.bash 
  541  cd RECYCLE_MODEL/recycle_model/
  542  cd STORE_model/
  543  git st
  544  cp  run_files/run_lsq_ponly_sbatch.bash run_files/hive/run_vpros_ponly_sbatch.bash run_files_1/
  545  rm run_files/.run_shgo_ponly_hive.bash.swp
  546  rm -rf VPRO/LSQ/
  547  mv run_files/copy_vpro.bash run_files/README.txt run_files_1/
  548  git checkout -- .
  549  git st
  550  git pull
  551  tmux ls
  552  tmux new -t monte
  553  tmux new -s monte
  554  cd RECYCLE_MODEL/recycle_model/
  555  cd STORE_model/
  556  ls run_files
  557  ls run_files/*
  558  ls run_files/hive/run_lsq_ponly_EXO.sbatch 
  559  cat run_files/hive/run_lsq_ponly_MIN.sbatch 
  560  cp run_files/hive/run_lsq_ponly_MIN.sbatch run_files/hive/run_monte_MIN.sbatch
  561  vim run_files/hive/run_monte_MIN.sbatch
  562  cat run_files/run_lsq_ponly_sbatch.bash 
  563  for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; mkdir /mnt/beegfs/home/
  564  for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do mkdir /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte_${m}; done
  565  for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do cp run_files/run_mo; done
  566  for m in MIN MIXOTROPH OVERFLOW ROS EXOENZYME; do perl -pe "s/MIN/${m}/g" run_files/hive/run_monte_MIN.sbatch > run_files/hive/run_monte_${m}.sbatch; done
  567  tail -n 1000 run_files/hive//run_monte_* | less
  568  cp run_files/hive/run_lsq_ponly_MIN.sbatch run_files/hive/run_monte_MIN.sbatch
  569  vim run_files/hive/run_monte_MIN.sbatch
  570  for m in MIXOTROPH OVERFLOW ROS EXOENZYME; do perl -pe "s/MIN/${m}/g" run_files/hive/run_monte_MIN.sbatch > run_files/hive/run_monte_${m}.sbatch; done
  571  tail -n 1000 run_files/hive//run_monte_* | less
  572  ls run_files/hive//run_monte_*
  573  cd slurm/
  574  ls ../run_files/hive//run_monte_*
  575  ../run_files/hive//run_monte_EXOENZYME.sbatch
  576  sbatch ../run_files/hive//run_monte_EXOENZYME.sbatch
  577  sbatch ../run_files/hive//run_monte_MIXOTROPH.sbatch 
  578  sbatch ../run_files/hive//run_monte_MIN.sbatch 
  579  sbatch ../run_files/hive//run_monte_ROS.sbatch 
  580  sbatch ../run_files/hive//run_monte_OVERFLOW.sbatch 
  581  squeue 
  582  ls
  583  ls -ltr
  584  cat ../run_files/hive//run_monte_OVERFLOW.sbatch 
  585  baseid=monte
  586  RDIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/recycle_model/STORE_model
  587  ODIR=/mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/${baseid}
  588  ls $ODIR 
  589  mkdir -p $ODIR/log
  590  sbatch ../run_files/hive//run_monte_OVERFLOW.sbatch 
  591  sbatch ../run_files/hive//run_monte_ROS.sbatch 
  592  sbatch ../run_files/hive//run_monte_MIN.sbatch 
  593  sbatch ../run_files/hive//run_monte_MIXOTROPH.sbatch 
  594  sbatch ../run_files/hive//run_monte_EXOENZYME.sbatch
  595  squeue 
  596  ls $ODIR
  597  tail /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/* | less
  598  tail /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/* 
  599  squeue  -u
  600  squeue  -u $USER | less
  601  cd ../
  602  ls ../../results/monte/out/
  603  ls ../../results/monte/out/|wc
  604  tmux attach
  605  squeue  -u $USER | less
  606  ls ../run_files/hive//run_monte_*
  607  cd RECYCLE_MODEL/recycle_model/STORE_model/
  608  ls ../../results/monte/out/|wc
  609  squeue  -u $USER | less
  610  sacct|less
  611  tail /mnt/beegfs/home/dsher/oweissber/RECYCLE_MODEL/results/monte/log/* | less
  612  ls ../../results/monte/out/|wc
  613  ls -ltr ../../results/monte/out/|less
  614  srun ls -ltr ../../results/monte/out/|less
  615  ls ../../results/monte/out/|wc
  616  bc 128661/3
  617  bc 
  618  ls ../../results/monte/out/|wc
  619  squeue  -u $USER 
  620  ls ../../results/monte/out/|wc
  621  squeue  -u $USER 
  622  ls -ltr ../../results/monte/log/|tail
  623  ../../results/monte/log/out_2788255_14_2788269.out 
  624  tail../../results/monte/log/out_2788255_14_2788269.out 
  625  less ../../results/monte/log/out_2788255_14_2788269.out 
  626  git st
  627  srun ./results_postprocess.py --dpath ../../results/monte/out/ --out results --run_id monte
  628   
  629  ls -ltrh results
  630  squeue  -u $USER 
  631  ls -ltrh results
  632  cd RECYCLE_MODEL/recycle_model/STORE_model/
  633  ls -ltrh results
  634  ls ../../results/monte/out/|wc
  635  ls -ltr ../../results/monte/out/|tail
  636  srun ls -ltr ../../results/monte/out/|tail
  637  srun ls -ltr ../../results/monte/out/monte_ROS_2_monte_vpro* | tail
  638  srun ls -ltr ../../results/monte/out/monte_ROS_36_monte_vpro* | tail
  639  srun ls -ltr ../../results/monte/out/monte_ROS_*_monte_vpro* | tail
  640  scancel 
  641  scancel -u $USER
  642  squeue  -u $USER 
  643  git st
  644  sbatch --wrap ' ./results_postprocess.py --dpath ../../results/monte/out/ --out results --run_id monte'
  645  squeue  -u $USER 
  646  ls -ltr
  647  less slurm*
  648  squeue  -u $USER 
  649  less slurm*
  650  less slurm-*
  651  cat slurm-*
  652  ls -ltr results
  653  ls -ltrh results
  654  srun cp -rp results/monte_* tmp/
  655  ls -ltrh results
  656  git st
  657  cat slurm-*
  658  squeue  -u $USER 
  659  sacct| less
  660  ls -ltrh results
  661  ls VPRO/LSQ/
  662  ls VPRO/LSQ/ROS/
  663  vim run_files/hive/run_monte_*
  664  git fetch
  665  mv run_files/hive/run_monte_* run_files_1/
  666  git pull
  667  vim run_files/hive/run_monte_*
  668  git diff
  669  git diff | grep -e '^\+'
  670  sbatch run_files/hive/run_monte_ROS.sbatch 
  671  sbatch run_files/hive/run_monte_OVERFLOW.sbatch 
  672  squeue 
  673  ls ../../results/monte2/out/
  674  ls ../../results/monte2/out/|wc
  675  ls ../../results/monte2/log/
  676  exit
  677  squeue 
  678  cd RECYCLE_MODEL/recycle_model/STORE_model/
  679  ls ../../results/monte2/log/
  680  ls ../../results/monte2/out/|wc
  681  tmux attach
  682  tmux attach
  683  sinfo
  684  tmux attach
  685  cd RECYCLE_MODEL/results/
  686  ls
  687  ls -d monte*
  688  ls -d monte*/out | less
  689  ls  monte*/out | less
  690  ls  monte3*/out/*  | less
  691  cd ../recycle_model/STORE_model/
  692  ls results
  693  squeue 
  694  sinfo
  695  tmux attach
  696  tmux new -t kmean
  697  tmux new kmean
  698  tmux new -a kmean
  699  tmux new -s kmean
  700  squeue 
  701  sacct
  702  squeue 
  703  sacct
  704  squeue 
  705  sacct
  706  squeue 
  707  sacct
  708  sacct | less
  709  squeue | less
  710  sacct | less
  711  tmux ls
  712  tmux attach
  713  tmux attach
  714  squeue 
  715  cd RECYCLE_MODEL/recycle_model/STORE_model/slurm/
  716  cat slurm-28499*
  717  squeue 
  718  cd ../
  719  srun ./results_postprocess.py --dpath ~/RECYCLE_MODEL/results/baseline/out/ --out results --run baseline_overflow_new
  720  tmux attach
  721  cd RECYCLE_MODEL/recycle_model/
  722  git st
  723  git checkout main
  724  srun git checkout main
  725  srun git pull
  726  cd STORE_model/
  727  ls
  728  ls results
  729  cd results
  730  mkdir old_results_25122023
  731  ls
  732  *.* old_results_25122023/
  733  mv *.* old_results_25122023/
  734  ls
  735  mv loss_analysis/ old_results_25122023/
  736  git st
  737  srun git st
  738  cd ..
  739  ls
  740  ls EXOENZYME/
  741  rm EXOENZYME/
  742  rmdir EXOENZYME/
  743  rmdir 
  744  rmdir MIXOTROPH/
  745  rmdir MIN/
  746  rmdir OVERFLOW/
  747  ls
  748  ls run_files
  749  ls  -ltr run_files
  750  srun rm -rf ~/RECYCLE_MODEL/results/*
  751  bg
  752  squeue 
  753  sacct
  754  ls ~/RECYCLE_MODEL/results/
  755  git log
  756  git pull
  757  git st
  758  less run_files/run_baseline.sbatch
  759  less run_files/run_sensitivity
  760  less run_files/run_sensitivity.sbatch 
  761  less run_files/run_baseline.sbatch
  762  ls slurm/
  763  srun rm -rf slurm/*
  764  jobs
  765  ls ~/RECYCLE_MODEL/results/
  766  ls run_files
  767  ls run_files/run_monte.bash 
  768  cat run_files/run_monte.bash 
  769  cat run_files/hive/run_monte_ROS.sbatch 
  770  vim run_files/run_sensitivity.sbatch 
  771  ls run_files/*99*
  772  cat run_files/run_sensitivity_pro99.bash
  773  vim run_files/run_sensitivity.sbatch 
  774  vim run_files/run_baseline.sbatch 
  775  cd slurm/
  776  source ../run_files/run_sensitivity.sbatch 
  777  source ../run_files/run_baseline.sbatch 
  778  sacct
  779  sacct | less -S 
  780  squeue 
  781  squeue -u $USER 
  782  cd RECYCLE_MODEL/results/sensitivity/out/
  783  ls
  784  tmux ls
  785  tmux attach 
  786  sacct
  787  sacct | less -S 
  788  squeue | less -S 
  789  sinfo 
  790  squeue | less -S 
  791  squeue -u $USER | grep PD | less -S
  792  squeue -u $USER | grep -v PD | less -S
  793  ls
  794  cd RECYCLE_MODEL/recycle_model/STORE_model/slurm/
  795  ls
  796  less * | less
  797  less * | grep MSE | less
  798  less * | grep MSE | sort -n | less
  799  sbatch --wrap 'less * | grep MSE | sort -n'
  800  less slurm-2950071.out 
  801  sacct | less -S 
  802  squeue -u $USER | grep -v PD 
  803  less slurm-2949476.out 
  804  less slurm-2950071.out 
  805  ls ~/RECYCLE_MODEL/results/
  806  ls ~/RECYCLE_MODEL/results/baseline/log/
  807  ls ~/RECYCLE_MODEL/results/baseline/
  808  ls ~/RECYCLE_MODEL/results/sensitivity/
  809  ls ~/RECYCLE_MODEL/results/sensitivity/out/
  810  ls ~/RECYCLE_MODEL/results/sensitivity/out/ |wc -l
  811   ~/RECYCLE_MODEL/results/sensitivity/out/ | less
  812  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | less
  813  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | perl -ne 'split("_"); print("$_[0] $_[6]\n")' | less
  814  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | perl -ne 'split("_"); print("$_[1] $_[6]\n")' | sort -u
  815  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | perl -ne 'split("_"); print("$_[1] $_[6]\n")' | uniq -c 
  816  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | perl -ne 'split("_"); print("$_[1] $_[6]\n")' | uniq -c | less
  817  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | perl -ne 'split("_"); print("$_[1] $_[6]\n")' | sort | uniq -c | less
  818  cat ../run_files/run_sensitivity.sbatch 
  819  cat ../run_files/run_baseline
  820  cat ../run_files/run_baseline.sbatch 
  821  sacct
  822  squeue -u $USER | grep -v PD 
  823  squeue -u $USER | grep PD 
  824  squeue -u $USER | grep PD |wc
  825  squeue -u $USER | grep -v PD |wc
  826  squeue -u $USER | grep PD 
  827  squeue -u $USER | grep -v PD |wc
  828  squeue -u $USER | grep PD |wc
  829  squeue -u $USER | grep -v PD 
  830  scancel 
  831  scancel -u $USER
  832  squeue -u $USER 
  833  sacct | less -S
  834  vim ../run_files/run_baseline.sbatch 
  835  source ../run_files/run_baseline.sbatch
  836  squeue -u $USER 
  837  vim ../run_files/run_sensitivity.sbatch 
  838  squeue -u $USER 
  839  ls  ~/RECYCLE_MODEL/results/baseline/
  840  ls  ~/RECYCLE_MODEL/results/baseline/out/
  841  ls  ~/RECYCLE_MODEL/results/baseline/out/|wc
  842  vim ../run_files/run_sensitivity.sbatch 
  843  sleep 1
  844  sleep 2
  845  vim ../run_files/run_sensitivity.sbatch 
  846  sinfo
  847  ls
  848  squeue -u $USER 
  849  rm -rf slurm-29*
  850  mv ~/RECYCLE_MODEL/results/sensitivity/ ~/RECYCLE_MODEL/results/sensitivity_problematic
  851  source ../run_files/run_sensitivity.sbatch 
  852  squeue -u $USER 
  853  squeue -u $USER | less
  854  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ | perl -ne 'split("_"); print("$_[1] $_[6]\n")' | sort | uniq -c | less
  855  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ 
  856  ls  ~/RECYCLE_MODEL/results/sensitivity
  857  ls
  858  less slurm-2950* | less
  859  tail slurm-2950* | less
  860  ls  ~/RECYCLE_MODEL/results/sensitivity/out/ 



#####################################################################################



  963  srun ./create_sbatch_for_multimodel.py --template run_files/hive/multi/template/run_monte_ponly_template --outdpath run_files/hive/multi/run_files/ponly/
 1014  srun ./create_sbatch_for_multimodel.py  --template run_files/hive/multi/template/run_monte_ponly.template --out run_files/hive/multi/run_files/ponly/
 1024  for f in ../run_files/hive/multi/run_files/ponly/*; do sbatch $f; done

 1048  mkdir -p ./results/multi/ponly;dp=/lustre1/home/dsher/oweissber/RECYCLE_MODEL/results/multi/monte_ponly_;for i in ~/RECYCLE_MODEL/results/multi/monte_ponly_*; do m=${i/${dp}/}; sbatch --wrap "./results_postprocess.py --dpath $i/out/ --outdpath ./results/multi/ponly/ --run_id monte_ponly_${m}" ; done



 1172  source ../run_files/hive/multi/run_find_vpro.bash



 1190  ./create_sbatch_for_multimodel.py --template run_files/hive/multi/template/run_monte_het.template  --outdpath run_files/hive/multi/run_files/het/
 1197  for i in  ../run_files/hive/multi/run_files/het/* ; do sbatch $i; done



 1234  srun ./create_sbatch_for_multimodel.py --template run_files/hive/multi/template/run_monte_ponly_extend.template --outdpath run_files/hive/multi/run_files/ponly_extend/ --vprodpath $PWD/VPRO/multi/run1/
 1243  for i in ../run_files/hive/multi/run_files/ponly_extend/*; do sbatch $i; done
 1252  srun ./create_sbatch_for_multimodel.py --template run_files/hive/multi/template/run_monte_ponly_extend.template --outdpath run_files/hive/multi/run_files/ponly_extend_single/ --vprodpath $PWD/VPRO/X0

 1256  for i in ../run_files/hive/multi/run_files/ponly_extend_single/*; do sbatch $i; done




##################

sbatch --wrap "./results_postprocess.py --dpath ~/RECYCLE_MODEL/results/multi/het/monte/out/  --outdpath ./results/multi/ --run_id monte_het_multi"
sbatch --wrap "./results_postprocess.py --dpath ~/RECYCLE_MODEL/results/multi/ponly_extend/monte_add_ponly_*/out/   --outdpath ./results/multi/ --run_id monte_extend_ponly"

