#!/usr/bin/env bash
echo "Not ready yet but it might work.  Continue? (y/n)"
case $A in
"n")
exit 1
;;
"y")
break
;;
*)
echo "y/n?"
;;

clear
adir=$HOME/ast07
if [ -d "$adir" ];
then
echo "$HOME/ast07 already exists.  What should I do?"
select action in "Backup $adir to $adir/bak and continue" "Write into $adir/ethan (breaks files!)" "Cancel"

#while true; 
do
case $action in
"Backup $adir to $adir/bak and continue")
mkdir -v $adir/bak
mv $adir/ $adir/bak
echo "(Ignore error ^)"
mv $(pwd) $adir/
mv $adir/AST-07/ast07/* $adir
rm -rf $adir/AST-07/ast07
#mv $(pwd)/* $adir/
#mkdir -v '$adir/AST-07'
#mv $adir/*.*py* $adir/AST-07
rm -rf $pwd
break
;;
"Write into $adir/ethan (breaks files!)")
echo "Creating $adir/ethan and moving files..."
mkdir -v $adir/ethan
mv $(pwd) $adir/ethan
break
;;
"Cancel")
echo "Cancel."
exit 1
;;
*)
echo "Select an option."
;;
esac
done

elif [ ! -d $adir ];
then
mv -v "$(pwd) $adir"
mv -v $adir/AST-07/ast07/* $adir/
rm -rfv $adir/AST-07/ast07
fi

