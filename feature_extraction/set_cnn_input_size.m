function net = set_cnn_input_size(net, im_size)

net.meta.normalization.imageSize(1:2) = round(im_size(1:2));
net.meta.normalization.averageImage = imresize(single(net.meta.normalization.averageImageOrig), net.meta.normalization.imageSize(1:2));

net.info = vl_simplenn_display(net);
end