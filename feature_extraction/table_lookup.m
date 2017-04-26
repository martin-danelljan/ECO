function out=table_lookup(im,table)

if isa(im,'uint8')
    im = int32(im);
end;

[im_height, im_width, num_im_chan, num_images] = size(im);

den = int32(8);
fac = int32(32);
offset = int32(1);

if num_im_chan == 3
    RR=im(:,:,1,:);GG=im(:,:,2,:);BB=im(:,:,3,:);
    index_im = offset + idivide(RR,den) + fac*idivide(GG,den) + fac*fac*idivide(BB,den);
    out = permute(reshape(table(index_im(:),:), im_height, im_width, num_images, size(table,2)), [1 2 4 3]);
else
    out = permute(reshape(table(im(:)+1,:), im_height, im_width, num_images, size(table,2)), [1 2 4 3]);
end;
