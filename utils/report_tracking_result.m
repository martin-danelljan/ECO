function seq = report_tracking_result(seq, result)

if strcmpi(seq.format, 'otb')
    seq.rect_position(seq.frame,:) = round([result.center_pos([2,1]) - (result.target_size([2,1]) - 1)/2, result.target_size([2,1])]);
elseif strcmpi(seq.format, 'vot')
    if seq.frame > 1
        bb_scale = 1;
        sz = result.target_size / bb_scale;
        tl = result.center_pos - (sz - 1)/2;
        br = result.center_pos + (sz - 1)/2;
        x1 = tl(2); y1 = tl(1);
        x2 = br(2); y2 = br(1);
        result_box = round(double([x1 y1 x2 y1 x2 y2 x1 y2]));
        if any(isnan(result_box) | isinf(result_box))
            error('Illegal values in the result.')
        end
%         if any(result_box < 0)
%             error('Negative values')
%         end
        seq.handle = seq.handle.report(seq.handle, result_box);
    end
else
    error('Uknown sequence format');
end