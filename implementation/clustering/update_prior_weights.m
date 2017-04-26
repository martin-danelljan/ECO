function [prior_weights, replace_ind] = update_prior_weights(prior_weights, sample_weights, latest_ind, frame, params)

% Update the training sample wights

learning_rate = params.learning_rate;

if frame == 1
    replace_ind = 1;
%     prior_weights = zeros(params.nSamples,1,'single');
    prior_weights(replace_ind) = 1;
else
    % Get which sample to replace.
    switch params.sample_replace_strategy
        case 'lowest_prior'
            [~, replace_ind] = min(prior_weights);
        case 'lowest_weight'
            [~,replace_ind] = min(sample_weights);
        case 'lowest_median_prior'
            median_prior = median(prior_weights);
            % add 1 to all the weights that has a prior weight larger
            % than the median, and then take the mininum
            [~,replace_ind] = min(sample_weights + (prior_weights > median_prior));
        case 'constant_tail'
            [~,I] = sort(prior_weights,'ascend');
            lt_index = I(1:params.lt_size);
            st_index = I(params.lt_size+1:end);
            
            %some things are best forgotten...
            dummy_w = sample_weights;
            dummy_w(st_index) = inf;
            [~, replace_ind] = min(dummy_w);
    end
    
    if frame == 2
        prior_weights(latest_ind) = 1 - learning_rate;
        prior_weights(replace_ind) = learning_rate;
    else
        % take the previous value and ensure a relative difference of
        % (1-learning_rate)
        prior_weights(replace_ind) = prior_weights(latest_ind) / (1 - learning_rate);
    end
    
    if params.lt_size > 0
        [~,I] = sort(prior_weights,'ascend');
        lt_index = I(1:params.lt_size);
        st_index = I(params.lt_size+1:end);
        
        minw = min(prior_weights(st_index));
        
        if minw ~= 0
            lt_mask = false(size(prior_weights));
            lt_mask(lt_index) = true;
            lt_mask = lt_mask & (prior_weights > 0);
            prior_weights(lt_mask) = minw * (1 - learning_rate);
        end;
    end;
    
    prior_weights = prior_weights/sum(prior_weights);
end