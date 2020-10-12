#! /bin/ksh

for FILE in *.fa;
do
if [[ -s $FILE ]] ; then
echo "$FILE" >> chr_genome_fa.txt
else
echo "$FILE" >> empty_file.txt
fi ;
done