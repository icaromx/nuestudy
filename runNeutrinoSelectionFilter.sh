RDIR='results'
mkdir ${RDIR}
#SPATH='/Users/ivan/Work/eLEE/nue1e0p_selection_v12/v17'

SAMPLES=(
	nue
	numu
	databnb
	dirt
	ext
	)

FILENAMES=(
	/Users/ivan/Work/eLEE/Samples/v14Samples/bnb_nue_cc0pinp/out.root
	/Users/ivan/Work/eLEE/Samples/v14Samples/bnb_nu_cc0pinp/out.root
	/Users/ivan/Work/eLEE/Samples/v14Samples/beam_on_cc0pinp/out.root
	/Users/ivan/Work/eLEE/Samples/v14Samples/bnb_dirt_cc0pinp/out.root
	/Users/ivan/Work/eLEE/Samples/v14Samples/beam_off_cc0pinp/out.root
)

#FILENAMES=(
#	/Users/ivan/Work/eLEE/Samples/nue1e0p_selection_v12/v17/prodgenie_bnb_intrinsic_nue_uboone_overlay_mcc9.1_run1_reco2.root
#	/Users/ivan/Work/eLEE/Samples/nue1e0p_selection_v12/v17/prodgenie_bnb_nu_uboone_overlay_mcc9.1_run1_reco2.root
#	/Users/ivan/Work/eLEE/Samples/nue1e0p_selection_v12/v17/data_bnb_optfilter_C1_5e19_goodruns_mcc9.1_reco2.root
#	/Users/ivan/Work/eLEE/Samples/nue1e0p_selection_v12/v17/prodgenie_bnb_dirt_overlay_run1_mcc9.1_v08_00_00_16_reco2.root
#	/Users/ivan/Work/eLEE/Samples/nue1e0p_selection_v12/v17/data_extbnb_mcc9.1_v08_00_00_16_run1_reco2.root
#)

for (( i = 0; i < 5; i++ )); do
	root -l -q 'NeutrinoSelectionFilter.C("'${FILENAMES[i]}'","'${SAMPLES[i]}'")'
	mv ${SAMPLES[i]}.root ${RDIR}/
done
