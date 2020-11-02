##############
# HOW MANY PARTS?
IDX1=$(shell seq 1 1 10)
NF=10
#####################
#STEP OF EACH PART
IDX2=25000
#######################
# HOW many layers
IDX3=$(shell seq 0 1 3)


warn:
	echo 'What would you like to do?'

all_1:
	make compile
	mkdir BIL_up_down_plots
	make make_dir
	make link
	make make_BOX
	make binder
	make spp_corr
	make deconvolve
	make all_2
all_2:
	make vsfg
	make tot



make_dir:
	$(foreach num, $(IDX1), \
		$(shell mkdir part_$(num))\
	)
	mkdir total

make link:
	$(foreach num, $(IDX1), \
	       	$(shell cd part_$(num)/.;\
		ln -s ../../interface_definition/interface.xyz;ln -s ../../interface_definition/grid_interface;ln -s ../../interface_definition/pos_rebuilt.xyz; ln -s ../../interface_definition/vel.xyz;\
		)\
	)

make_BOX:
	$(foreach num, $(IDX1), \
		$(shell cd part_$(num)/.;\
			sed 's/pippo1/$(shell echo "(($(num)-1)*$(IDX2)+1)" | bc -l)/' ../BOXtemp | sed 's/pippo2/$(shell echo "(($(num))*$(IDX2))" | bc -l)/'  > BOXDATA;\
		)\
	)			


binder:
	$(foreach num, $(IDX1), \
	       	$(shell cd part_$(num)/.; ../create_binder_istart.x >> out.log;\
		)\
	)

spp_corr:
	$(foreach num, $(IDX1), \
	$(foreach num2, $(IDX3), \
	       	$(shell cd part_$(num)/.; ../ssp_sfg_selfcorr_istart.x binder_dL$(num2).xyz corr_dL$(num2).dat >> out.log ;\
		)\
	)\
	)

vsfg:
	$(foreach num, $(IDX1), \
	$(foreach num2, $(IDX3), \
	       	$(shell cd part_$(num)/.; ../vsfg_signal_morita.x corr_dL$(num2).dat ssp_dL$(num2).dat >> out.log;\
		)\
		)\
	)


deconvolve:
	$(foreach num, $(IDX1), \
		$(shell cd part_$(num)/.; ../deconvolve.x >> out.log; ../ssp_sfg_selfcorr_istart.x binder_down.xyz corr_down.dat >> out.log; ../vsfg_signal_morita.x corr_down.dat ssp_down.dat >> out.log; ../ssp_sfg_selfcorr_istart.x binder_up.xyz corr_up.dat >> out.log; ../vsfg_signal_morita.x corr_up.dat ssp_up.dat >> out.log;\
		)\
	)

	

tot:
	$(foreach num, $(IDX1), \
	$(foreach num2, $(IDX3), \
		 $(shell cd part_$(num)/.; cp ssp_dL$(num2).dat ../total/ssp_dL$(num2)_$(num) ; cp ssp_down.dat ../total/ssp_down_$(num); cp ssp_up.dat ../total/ssp_up_$(num)\
		)\
		)\
		 $(shell cd part_$(num)/.;  ../sum_spect.x ssp_dL0.dat ssp_dL1.dat ssp_BIL.dat;../sum_spect.x ssp_BIL.dat ssp_dL2.dat TOT.dat; \
		)\
	)

sum:
	$(foreach num2, $(IDX3), \
		 $(shell cd total/.; ../sum_part.x $(NF) ssp_dL$(num2)_  ssp_dL$(num2)_tot.dat;\
		)\
	)
	cd total/.; ../sum_spect.x ssp_dL0_tot.dat ssp_dL1_tot.dat ssp_BIL_tot.dat; cd ..
	cd total/.; ../sum_spect.x ssp_BIL_tot.dat ssp_dL2_tot.dat Imchi2.dat; cd ..


sum2:
	$(foreach num2, $(IDX3), \
		 $(shell cd total2/.; ../sum2_part.x 7 3  ssp_dL$(num2)_  ssp_dL$(num2)_tot.dat;\
		)\
	)
	cd total2/.; ../sum_spect.x ssp_dL0_tot.dat ssp_dL1_tot.dat ssp_BIL_tot.dat; cd ..
	cd total2/.; ../sum_spect.x ssp_BIL_tot.dat ssp_dL2_tot.dat Imchi2.dat; cd ..

kill:
	rm -r part_*
	rm -r total

compile:
	gfortran -o ssp_sfg_selfcorr_istart.x ssp_sfg_selfcorr_istart.f90 
	gfortran -o create_binder_istart.x create_binder_istart.f90 
	gfortran -o vsfg_signal_morita.x vsfg_signal_morita.f
	gfortran -o deconvolve.x deconvolve.f90

updownPlot:
	$(foreach num, $(IDX1), \
		$(shell cd part_$(num)/.; python3.8 ../SFG_plotter.py ssp_up.dat ssp_down.dat ssp_dL1.dat L1_upDownTot_part_$(num).png ; cp L1_upDownTot_part_$(num).png ../BIL_up_down_plots; cd ..;\
		)\
	)

