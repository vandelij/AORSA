import os, sys
import time

with open('shot_number.txt') as f:
    lines = f.readlines()

shotnum = lines[0] #'testing' # lines[0]
time_stamp = round(time.time())


#these shenanigans relate to vscode not having the working directory as the directory of the file it runs
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

host = 'perlmutter-p1.nersc.gov'
username = 'vandelij'
parent_direc = '~/aorsa/test/diiid-aorsa-hires_copy'# '/global/homes/j/jwright/perlmutter-builds/aorsa/examples/DIIID-helicon'#'~/AORSA/DIIID-helicon/'  #'~/AORSA/DIIID-helicon/'
file_to_grab = 'aorsa2d.ps' #'aorsa2d_input_og_dont_tuch.in'
file_Efield_2D = 'Efield_2D.vtk'
target_direc = '~/Desktop/HFS_HHFW_antenna/AORSA_MCGO_Scripts/AORSA/shots/'
file_to_save = f'{shotnum}/aorsa{time_stamp}.ps'#f'{shotnum}/aorsa{time_stamp}.ps' 
# os.system(f'scp {username}@{host}:{rwdir}/cql3d.nc shots/{shotNum}/cql3d.nc')
# os.system(f'scp {username}@{host}:{rwdir}/cql3d_krf001.nc shots/{shotNum}/cql3d_krf001.nc')

os.system(f'scp {username}@{host}:{parent_direc}/{file_to_grab} {target_direc}{file_to_save}')
os.system(f'scp {username}@{host}:{parent_direc}/{file_Efield_2D} {target_direc}/{shotnum}/Efield_2D_new.vtk')

# host = 'eofe7.mit.edu'
# username = 'vandelij'
# parent_direc = '/home/vandelij/HFW_174658/case_beampwr_10_general_electron_general_ion_colmodl0_both_rf_copy_to_fix_power_issue/'# '/global/homes/j/jwright/perlmutter-builds/aorsa/examples/DIIID-helicon'#'~/AORSA/DIIID-helicon/'  #'~/AORSA/DIIID-helicon/'
# file_to_grab = 'genray.in' #'aorsa2d_input_og_dont_tuch.in'
# file_Efield_2D = 'Efield_2D.vtk'
# target_direc = '~/Desktop/HFS_HHFW_antenna/AORSA_MCGO_Scripts/AORSA/shots/'
# file_to_save = f'{shotnum}/genray.in'#f'{shotnum}/aorsa{time_stamp}.ps' 
# # os.system(f'scp {username}@{host}:{rwdir}/cql3d.nc shots/{shotNum}/cql3d.nc')
# # os.system(f'scp {username}@{host}:{rwdir}/cql3d_krf001.nc shots/{shotNum}/cql3d_krf001.nc')

# os.system(f'scp {username}@{host}:{parent_direc}{file_to_grab} {target_direc}{file_to_save}')




#os.system(f'scp {username}@{host}:~/AORSA/HFW_174658/aorsa2d.ps ~/Desktop/HFS_HHFW_antenna/AORSA_MCGO_Scripts/AORSA/shots/{shotnum}/aorsa{time_stamp}.ps')
#os.system(f'scp {username}@{host}:~/AORSA/DIIID-helicon/DIII_NB_FW_0th.max.nc ~/Desktop/HFS_HHFW_antenna/AORSA_MCGO_Scripts/AORSA/DIII_NB_FW_0th.max.nc')
#os.system(f'scp {username}@{host}:~/AORSA/DIIID-helicon/{file_to_grab} ~/Desktop/HFS_HHFW_antenna/AORSA_MCGO_Scripts/AORSA/')