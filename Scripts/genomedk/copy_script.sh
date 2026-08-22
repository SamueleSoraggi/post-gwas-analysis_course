FOLDERS=`find ../ -maxdepth 1 -type d ! -name test_samuele ! -name input ! -path "../"`

for F in $FOLDERS
do
	echo $F
	cp gwas-minimal.sif $F/
	cp start $F/
	ln -sf /faststorage/project/post-gwas-phd-course/data/input $F/input
        ln -sf /faststorage/project/post-gwas-phd-course/data/reference_data $F/reference_data
	cp *.ipynb $F/
	cp -a scripts $F/
	cp -a figures $F/
	chmod -R 777 $F
	echo "done $F"
	
done
