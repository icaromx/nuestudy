RDIR='results'
#mkdir ${RDIR}
#SPATH='/Users/ivan/Work/eLEE/nue1e0p_selection_v12/v17'



SAMPLES=(
	nue
	numu
	databnb
	dirt
	ext
	)

iter=5
input=$2
if [ "$1" == 'test' ]; then
#        iter=1
				input="/Users/ivan/Work/eLEE/FileLists/v12.list"
fi

i=0
while IFS= read -r line
do
	root -l -q 'NeutrinoSelectionFilter.C("'$line'","'${SAMPLES[i]}'","'$1'")'
	((i+=1))
	if [ "$iter" == '1' ]; then
	        break
	fi
done < "$input"


#for (( i = 0; i < $iter; i++ )); do
#	root -l -q 'NeutrinoSelectionFilter.C("'${FILENAMES[i]}'","'${SAMPLES[i]}'","'$1'")'
#	mv ${SAMPLES[i]}.root ${RDIR}/
#done
