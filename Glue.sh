#!/bin/bash
#args: filelist targetname
i=0
list_aux=`echo "  "`
list=$list_aux
for file in `cat $1`
do
  list=`echo -n "$list_aux  $file"` 
  list_aux=$list
  let i+=1
done
echo "Concatenating $i files..."
cat $list > $2
echo "Done."
