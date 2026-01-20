# --- The sending script for the Department of Geophysics cluster ---
# --- Target machine, user and directory ---
machine="nfsy6"
user="kihoulou"
dir="Titan_clathrates" 

# --- Launching script ---
cp run_MilleFEuiIle.sh /$machine/$user/$dir

# --- The parameter file ---
cp m_parameters_Titan.py /$machine/kihoulou/$dir

# --- MilleFEuiIle files ---
cp m_boundary_conditions.py /$machine/kihoulou/$dir
cp m_check.py /$machine/kihoulou/$dir
cp m_constants.py /$machine/kihoulou/$dir
cp m_elements.py /$machine/kihoulou/$dir
cp m_equations.py /$machine/kihoulou/$dir
cp m_filenames.py /$machine/kihoulou/$dir
cp m_interpolation.py /$machine/kihoulou/$dir
cp m_incompatibility.py /$machine/kihoulou/$dir
cp m_material_properties.py /$machine/kihoulou/$dir
cp m_melting.py /$machine/kihoulou/$dir
cp m_mesh.py /$machine/kihoulou/$dir
cp m_postproc.py /$machine/kihoulou/$dir
cp m_rheology.py /$machine/kihoulou/$dir
cp m_timestep.py /$machine/kihoulou/$dir
cp m_tracers.py /$machine/kihoulou/$dir
cp main.py /$machine/kihoulou/$dir
