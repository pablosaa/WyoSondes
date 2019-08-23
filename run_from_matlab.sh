#!/bin/bash

setenv USER pga082
setenv HOME /home/pga082
setenv PATH ${PATH}:/home/pga082/GFI/Repos/WyoSondes

/usr/local/bin/matlab -nodesktop -nodisplay -r "cd /home/pga082/GFI/Repos/WyoSondes; [data, metvar]= RASOBS_DOWNLOAD_DATA_RAW('04339',[2014:2019],[1,2,3,11,12],[1:31],[0 12], 'matfile',true,'waiting',5);exit" > pippo&

# [data, metvar]=RASOBS_DOWNLOAD_DATA_RAW" > pippo&
