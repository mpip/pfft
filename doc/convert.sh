files=*.tex
for i in $files; do
   # cat $i | sed 's/\\code{\([^}]*\)}/{\\small \1}/g' | 
    echo $i
    cat preample.tex $i | sed \
        -e 's;\hdots;\dots;g' \
        -e 's/\\code{\([^}]*\)}/\\verb+\1+/g' \
        | pandoc -f latex -t rst  \
    > ${i//.tex}.rst
done

