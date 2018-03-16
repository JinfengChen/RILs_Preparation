echo "summary PCR screen resuls"
python SV_screen_PCR_results.py

echo "Need PCR, missing data and SV numbers"
paste ../../Prepare0_mPing_SV_screen/bin/HEMP8_HEMP9.check.list HEMP8_HEMP9.PCR.results.rils.txt | perl -e 'while(<>){if ($_=~/RIL\d+,NA/ or $_=~/not_done/){print $_}}' | wc -l
paste ../../Prepare0_mPing_SV_screen/bin/HEMP8_HEMP9.check.list HEMP8_HEMP9.PCR.results.rils.txt | perl -e 'while(<>){if ($_=~/RIL\d+,NA/ or $_=~/not_done/){print $_}}' | cut -f2,3 | grep "na" | wc -
paste ../../Prepare0_mPing_SV_screen/bin/HEMP8_HEMP9.check.list HEMP8_HEMP9.PCR.results.rils.txt | perl -e 'while(<>){if ($_=~/RIL\d+,NA/ or $_=~/not_done/){print $_}}' | cut -f2,3 | grep "SV" | wc -ll
