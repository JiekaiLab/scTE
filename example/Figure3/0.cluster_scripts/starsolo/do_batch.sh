

for f in  ../fqs/*.p1.fq.gz
do
    root=`basename $f`
    path=`dirname $f`
    
    bf=`echo $root | sed -r 's#.p1.fq.gz##g'`
    p2=`echo $f | sed 's#.p1.fq.gz#.p2.fq.gz#g'`
    tt=`echo ss.$bf.Aligned.sortedByCoord.out.bam` # outfile
    if [ ! -f $tt ] # Check not already done
    then
        echo STARsolo $tt
        qsub -N solo.$bf -v p1=$f,p2=$p2,out=$bf. starsolo.sh
        sleep 2
    fi
done

