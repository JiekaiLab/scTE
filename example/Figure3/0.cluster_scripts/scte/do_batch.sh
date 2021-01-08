

for f in  ../starsolo*/*.bam
do
    root=`basename $f`
    path=`dirname $f`
    
    bf=`echo $root | sed -r 's#.Aligned.sortedByCoord.out.bam##g' | sed 's#.bam##g'`
    tt=`echo $bf.csv.gz ` # outfile
    if [ ! -f $tt ] # Check not already done
    then
        echo scTE $tt
        qsub -N scte.$bf -v in=$f,out=$bf scte.sh
        sleep 1
    fi
done

