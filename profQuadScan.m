data = E200_load_data('/Volumes/PWFA_4big/nas/nas-li20-pm01/E200/2013/20131108/E200_11469/E200_11469_scan_info.mat')

% Get an image
CEGAIN = data.raw.images.CEGAIN;

step6uid = E200_api_getUID(data.raw.scalars.step_num,6);
wantedUID = intersect(step6uid,CEGAIN.UID);
img = E200_load_images(CEGAIN,wantedUID(17));

% Get yvec for energy
[ysize,xsize] = size(img);
yvec = 1:ysize;
res = 10.3934;
e_axis = E200_cam_E_cal(data,yvec,res);
e_axis = fliplr(e_axis);
	
	
res = E200_api_getdat(CEGAIN,wantedUID(17),'RESOLUTION');
imagesc(img{1})

transfer.res = res;
transfer.e_axis = e_axis;

% Send to python

% ====================================
% Save for python loading
% ====================================
curpath = pwd();
savefile = fullfile(curpath,'tempfiles','forpython.mat');
save(savefile,'img','transfer','-v7');

