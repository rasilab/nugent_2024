for prefix in `cat sample_annotations.csv  | cut -f7 -d,`
do 
    files=`ls -1 ../*$prefix"_"*.fastq.gz`
    for file in $files
    do 
        ln -s $file .
    done
done
