for i in $(find . -type d) ; do 
    echo $i ; 
    ( find $i -type f | wc -l ) ; 
done