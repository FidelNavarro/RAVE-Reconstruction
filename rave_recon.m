rawfile = fullfile('/project','PRISM', 'fnav0118', 'Code_Fidel',  'PRJ-PRISM', 'RAVE_Tests', 'Volunteer2', 'TWIX', 'meas_MID00096_FID182203_RAVE.dat' );
output_path = fullfile('/project', 'PRISM', 'fnav0118', 'Code_Fidel', 'PRJ-PRISM','TWIX_outputs_007');
mode_file = fullfile('/project', 'PRISM', 'fnav0118', 'Code_Fidel', 'PRJ-PRISM', 'Code', 'yarra-dev-yarramodules-grasp-basic-6b525c658d42', 'modes', 'GRASP_basic.mode');
%spokes = [13, 21, 34, 55, 84, 89, 100] %this can  be use to reconstruct %(1)
%for several value of spokes. Special care with the mode file might be
%needed
%Regular function line to reconstruct using the mode_file instructions
yarra_GRASP_basic_FN_3(rawfile,output_path,mode_file); %(2)
%Function to reconstruct for different number of spokes using the "spokes"
%list
%yarra_GRASP_basic_FN_3(rawfile,output_path,mode_file, spokes); %(1)
