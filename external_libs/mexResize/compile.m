
% Compiles the mexResize function.

% Try to use the included mex files first. If they do not work, you can try
% to compile it yourself with this script. If the bellow mex command does
% not work, try to link to your own OpenCV installation.

if ispc
    mex -lopencv_core242 -lopencv_imgproc242 -L./ -I./ mexResize.cpp MxArray.cpp
else
    mex -lopencv_core -lopencv_imgproc -L./ -I./ mexResize.cpp MxArray.cpp
end