qsub -sync y -l h_vmem=10000M run_main.sh snets.params -g -ll 0 -lf log.txt
# qsub -sync yes -l h_vmem=3000M run_main.sh snets.params -g -ll 0 -lf log.txt
# qsub -sync yes -l h_vmem=3000M -t 1-10 run_step.sh ExecFilterPairs
# qsub -sync yes -l h_vmem=3000M run_main.sh snets.params -z -i filterpairs
