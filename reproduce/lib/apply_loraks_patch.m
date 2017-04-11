%Apply patch to LORAKS.m to 
%generate the file LORAKS_mod_GO.m
%that incorporates per iteration 
%error computations.
%This is necessary to reproduce 
%the plots in Figure 9 in the GIRAF paper.
%Requires UNIX 'patch' utility.
%Ensure LORAKS.m is in the MATLAB path.
%LORAKS_mod_GO.m will be added to the
%folder containing LORAKS.m

if ~exist('LORAKS.m','file')
    error('LORAKS.m not found');
end
originalfile = which('LORAKS.m');
[loraks_path,~,~] = fileparts(originalfile);
updatedfile = [loraks_path,'/LORAKS_mod_GO.m'];
patchfile = which('LORAKS.patch');

cmd = sprintf('patch %s -i %s -o %s',originalfile,patchfile,updatedfile);
[status,cmdout] = system(cmd,'-echo');
if status == 0 && exist('LORAKS_mod_GO.m','file')
    fprintf('sucessfully generated LORAKS_mod_GO.m\n');
else
    fprintf('failed to generate LORAKS_mod_GO.m\n');
end