#!/bin/bash

for start1 in `seq 3 13`
do

end1=14

    for num in `seq $start1 $end1`
        do

        cat Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_150mM_upreg_intnums_${num}.txt\
        >> Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_150mM_upreg_intnums_sep${start1}-${end1}.txt
    done

weblogo -f Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_150mM_upreg_intnums_sep${start1}-${end1}.txt\
 -D clustal -o intnum_150mM_${start1}-${end1}.png -F png -A dna -a ATGC -c classic -s large -n 8 --errorbars YES -i -2 --show-xaxis NO\
 -t "splice site conservation in genes with intron number ${start1}-${end1}, 150 mM" --title-fontsize 6

done

for end2 in `seq 2 13`
do

start2=1

    for num in `seq $start2 $end2`
        do

        cat Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_150mM_upreg_intnums_${num}.txt\
        >> Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_150mM_upreg_intnums_sep${start2}-${end2}.txt
    done

weblogo -f Arabidopsis_thaliana.TAIR10.34_isoform_SS_8_salt_150mM_upreg_intnums_sep${start2}-${end2}.txt\
 -D clustal -o intnum_150mM_${start2}-${end2}.png -F png -A dna -a ATGC -c classic -s large -n 8 --errorbars YES -i -2 --show-xaxis NO\
 -t "splice site conservation in genes with intron number ${start2}-${end2}, 150 mM" --title-fontsize 6

done
# other loop for the pair, then montage

for num in `seq 2 13`
do

num2=$((num + 1))

if [ "$num" != 13 ]; then
montage -mode concatenate -tile x1  intnum_150mM_1-${num}.png intnum_150mM_${num2}-14.png intnum_150mM_pair${num}.png
else
montage -mode concatenate -tile x1  intnum_150mM_1-${num}.png intnum_150mM_14.png intnum_150mM_pair${num}.png
fi

done
