P=psl
A=($P.ceftazidime $P.ceftriaxone $P.cefepime $P.tobramycin $P.gentamicin $P.cefazolin $P.ampicillin)

for i in ${A[@]}; do
    head -10 $1/consolidated/$i > $1/10/$i
    head -50 $1/consolidated/$i > $1/50/$i
    head -100 $1/consolidated/$i > $1/100/$i
    head -500 $1/consolidated/$i > $1/500/$i
    head -1000 $1/consolidated/$i > $1/1000/$i
    head -5000 $1/consolidated/$i > $1/5000/$i
    head -10000 $1/consolidated/$i > $1/10000/$i
done
