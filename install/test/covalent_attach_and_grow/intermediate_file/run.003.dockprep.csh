

  set mountdir = `pwd`

set dock = /nfs/home/tbalius/zzz.programs/dock6.6/dock6/bin

# ligand.sph is the whole ligand
cat << EOF > showbox.in 
Y                        
8.0                      
ligand.sph   
1                        
rec.box.pdb   
EOF

  $dock/showbox < showbox.in

# if (${name_prefix} != "rec") then 
#    ln -s ${mountdir}/rec/rec.box.pdb .
# endif

 set GRID = "$dock/grid"
ln -s /nfs/home/tbalius/zzz.programs/dock6.6/dock6/parameters/vdw_AMBER_parm99.defn .

echo "make grid"

#box_file                       ${name_prefix}.box.pdb
cat << EOF > ! grid.rec.in
compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             9
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  rec_charged_mod.mol2
box_file                       rec.box.pdb
vdw_definition_file            vdw_AMBER_parm99.defn
score_grid_prefix              rec.grid
EOF

${GRID} -i grid.rec.in -o grid.rec.out

cat << EOF > ! grid.rec.6_12.in
compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  rec_charged_mod.mol2
box_file                       rec.box.pdb
vdw_definition_file            vdw_AMBER_parm99.defn
score_grid_prefix              rec.6_12.grid
EOF

${GRID} -i grid.rec.6_12.in -o grid.rec.6_12.out


