function resizeddft  = resizeDFT(inputdft, desiredLen)

% Resizes a one-dimensional DFT to the desired length.

len = length(inputdft);
minsz = min(len, desiredLen);

scaling = desiredLen/len;

if size(inputdft, 1) > 1
    newSize = [desiredLen 1];
else
    newSize = [1 desiredLen];
end

resizeddft = complex(zeros(newSize, 'single'));

mids = ceil(minsz/2);
mide = floor((minsz-1)/2) - 1;

resizeddft(1:mids) = scaling * inputdft(1:mids);
resizeddft(end - mide:end) = scaling * inputdft(end - mide:end);
end