dirs=$(ls "./")
for dir in $dirs
do
	echo $dir
	if [ $dir != "clear.sh" ]
	then
		files=$(ls "$dir/")
		for file in $files
		do
			rm  "./$dir/$file"
		done

	fi
done
