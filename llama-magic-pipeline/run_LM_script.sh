#!/bin/bash
#SBATCH --time=5-00:00:00	
#SBATCH --mem=200G

p1=$1
p2=$2
p3=$3
p4=$4
p5=$5
p6=$6
p7=$7
p8=$8
p9=$9
p10=${10}
p11=${11}
p12=${12}
p13=${13}
p14=${14}
p15=${15}
p16=${16}
p17=${17}
p18=${18}
p19=${19}
p20=${20}

module load tandem
perl run_LM_script_pipeline.pl $p1 $p2 $p3 $p4 $p5 $p6 $p7 $p8 $p9 $p10 $p11 $p12 $p13 $p14 $p15 $p16 $p17 $p18 $p19 $p20