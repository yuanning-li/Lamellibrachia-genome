for FILENAME in *.tsv
do
echo $FILENAME >> $FILENAME.count
grep PF00452 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF00653 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF00319 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF00656 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF00531 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF01335 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF05729 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF00931 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
grep PF00020 $FILENAME | awk '{print $1}' | sort | uniq  | wc -l >> $FILENAME.count
done
