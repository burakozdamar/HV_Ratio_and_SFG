#########################
# BOX DATA
#########################

A= 13.386
B= 13.286
C= 85.000

#---------------------------
############################

INPUTPOS= pos.xyz


OUTPUT= pos_recentered.xyz

############################
#---------------------------
#############################

warn:
	echo 'What would you like to do?'


#############################
# COMPILE:
##############################
start: compile 1
compile:
	gfortran  ../source/Recenter_CM_periodic.f90 -o ../source/Recenter_CM_periodic.x
	gfortran  ../source/interface.f90 -o ../source/interface.x
	gfortran  ../source/rebuild_box.f90 -o ../source/rebuild_box.x
	gfortran  ../source/rebuild_box_solid.f90 -o ../source/rebuild_box_solid.x
	gfortran -o ../source/density_vs_r.x ../source/density_vs_r.f90
	gfortran -o ../source/up_down_water_layers.x ../source/up_down_water_layers.f90
	gfortran -o ../source/down_HBDONOR_network.x ../source/down_HBDONOR_network.f90
	gfortran -o ../source/down_hbond.x ../source/down_hbond.f90
	gfortran -o ../source/make_2D_graph_for_HB_network.x ../source/make_2D_graph_for_HB_network.f90
	gfortran -o ../source/PEG2Water-HB.x ../source/PEG2Water-HB.f90
	gfortran -o ../source/mapWC-Gaussian.x ../source/mapWC-Gaussian_istart.f90
###########################
# RUN:
#############################

1:
	make recenter
	make rebuild_box
	make rebuild_box_solid
2:
	make interface

1per10:
	./1per10pos.x

extract:
	./extractWater.x

recenter:
	$(shell	../source/Recenter_CM_periodic.x $(INPUTPOS) $(OUTPUT) $(A) $(B) $(C),  \
        )
rebuild_water:
	./../source/rebuild_water.x

rebuild_mol: 
	./../source/rebuild_molecules.x  

rebuild_box:
	 ./../source/rebuild_box.x  
rebuild_box_solid:
	 ./../source/rebuild_box_solid.x  
rebuild_box_KCl:
	./../source/rebuild_box_KCl.x 
      

interface:
	./../source/interface.x 

x=$(shell head -1 pos_rebuilt.xyz | xargs)
k=$(shell bc <<< $(x)\+2)

water_check:
	$(shell tail -n $k pos_rebuilt.xyz > t.xyz)
	$(shell head -n $k pos_rebuilt.xyz > h.xyz)
	@echo "\n"
	@echo "----First Step----"
	@python water_check.py h.xyz
	@echo "-----------------"
	@echo "\n"
	@echo "----Last Step----"
	@python water_check.py t.xyz
	@rm -fr t.xyz h.xyz
	@echo "\n"
