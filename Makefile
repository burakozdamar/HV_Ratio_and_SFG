#Cut the long traj into how many pieces?
IDX1=$(shell seq 1 1 10)
#Step of each part(to chang according to your traj)
IDX2=100
#Define layer number
NLAY=4


folders:
	make dir 
	make link
density:
	$(shell cd density_tot/.; ../source/density_vs_r.x; mv down_density_vs_r down_density_part_tot;)
BOX:
	$(foreach num, $(IDX1), $(shell cd part_$(num)/.; sed 's/dylan1/$(shell echo "(($(num)-1)*$(IDX2)+1)" | bc -l)/' ../BOXtemp | sed 's/dylan2/$(shell echo "(($(num))*$(IDX2))" | bc -l)/' > BOXDATA))
Hbonds:
	make layer
	make surfaceHb

dir:
	$(shell mkdir data density_tot)
	$(foreach num, $(IDX1), $(shell mkdir part_$(num)))


link:
	$(shell cd density_tot/.; ln -s ../interface_definition/interface.xyz; ln -s ../interface_definition/grid_interface; ln -s ../interface_definition/pos_rebuilt.xyz; ln -s ../interface_definition/pos_rebuilt_solid.xyz; ln -s ../interface_definition/BOXDATA) 
	$(foreach num, $(IDX1), $(shell cd part_$(num)/.; ln -s ../interface_definition/interface.xyz;ln -s ../interface_definition/grid_interface;ln -s ../interface_definition/pos_rebuilt.xyz ;ln -s ../interface_definition/pos_rebuilt_solid.xyz))

layer:
	$(foreach num, $(IDX1), $(shell cd part_$(num)/.; ../source/up_down_water_layers.x 4 ; ../source/down_hbond.x ; cp down_hbond ../data/restricted_down_hbond_part_$(num); mv up_down_water_layers.xyz up_down_restricted_part_$(num).xyz))


surfaceHb:
	$(foreach num, $(IDX1), $(shell cd part_$(num); ../source/PEG2Water-HB.x; cp PEG2Water.dat ../data/PEG2Water_part_$(num) ))

part_map:
	$(foreach num, $(IDX1), $(shell cd part_$(num)/.; ../source/mapWC-Gaussian.x ; python ../source/mapWlPlotter.py C-3-Gaussian-Map.dat; python ../source/mapWlPlotter.py O-3-Gaussian-Map.dat )) 
full_map:
	$(shell cd interface_definition; ../source/mapWC-Gaussian.x)
clean:
	rm -fr source/*.x 
	rm -fr part_*/ 
	rm -fr data/
	rm -fr density_tot/
	rm -fr interface_definition/pos_recentered.xyz
	rm -fr interface_definition/pos_rebuilt.xyz
	rm -fr interface_definition/pos_rebuilt_solid.xyz
	rm -fr interface_definition/*.dat
	rm -fr interface_definition/interface.xyz
