% data = E200_load_data('/Volumes/PWFA_4big/nas/nas-li20-pm01/E200/2013/20131108/E200_11469/E200_11469_scan_info.mat')

% Get an image
CEGAIN = data.raw.images.CEGAIN;

% Build shots I *don't* want
% Shot 1 of step 6
uids=E200_api_getUID(data.raw.scalars.step_num,6);
badUID = uids(1);
% badUID = 0;

% Subtract bad UIDs
wantedUIDs=CEGAIN.UID;
wantedUIDs=setdiff(wantedUIDs,badUID);

data = averageimages(data,'CEGAIN',wantedUIDs)

img=data.processed.images.CEGAIN.dat;

% Get yvec for energy
[ysize,xsize] = size(img{1});
yvec = 1:ysize;
res = 10.3934;

curpath=pwd();
addpath(fullfile(curpath,'E200_Cam_Energy_Cal'));

e_axis = E200_cam_E_cal(data,yvec,res);
e_axis = fliplr(e_axis);
	
	
res = E200_api_getdat(CEGAIN,wantedUIDs(17),'RESOLUTION');
imagesc(img{1})

transfer.res = res;
transfer.e_axis = e_axis;

% Send to python

% ====================================
% Save for python loading
% ====================================
curpath = pwd();
savefile = fullfile(curpath,'tempfiles','forpython.mat');
save(savefile,'data','transfer');

