for FileName in *.dat
	do
	head -n 3 $FileName | grep "^CLASSIFICATION:" | cut -d 
	"" -f 2,3
	done
