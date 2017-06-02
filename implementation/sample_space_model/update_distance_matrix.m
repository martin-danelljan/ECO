function [distance_matrix, gram_matrix]= update_distance_matrix(distance_matrix, gram_matrix, gram_vector, new_sample_norm, id1, id2, w1, w2)
% Updates the distance matrix

% Normalise the weights so that they sum to one
alpha1 = w1/(w1+w2);
alpha2 = 1 - alpha1;


if id2 < 0
    norm_id1 = gram_matrix(id1, id1);
    
    % Update the gram matrix
    if alpha1 == 0
        % The new sample replaces an existing sample.
        gram_matrix(:,id1) = gram_vector;
        gram_matrix(id1,:) = gram_matrix(:,id1);
        gram_matrix(id1, id1) = new_sample_norm;
    elseif alpha2 == 0
        % The new sample is discared 
    else
        % The new sample is merge with an existing sample
        gram_matrix(:,id1) = alpha1*gram_matrix(:,id1) + alpha2*gram_vector;
        gram_matrix(id1,:) = gram_matrix(:,id1);
        gram_matrix(id1, id1) = alpha1^2*norm_id1 + alpha2^2*new_sample_norm + 2*alpha1*alpha2*gram_vector(id1);
    end
    
    % Update distance matrix
    distance_matrix(:,id1) = max(gram_matrix(id1, id1) + diag(gram_matrix) - 2*gram_matrix(:,id1),0);
    distance_matrix(id1,:) = distance_matrix(:,id1) ;
    distance_matrix(id1,id1) = inf;
else
    if alpha1 == 0 || alpha2 == 0
        error('Error!');
    end
    
    % Two existing samples are merged and the new sample fills the empty
    % slot
    norm_id1 = gram_matrix(id1, id1);
    norm_id2 = gram_matrix(id2, id2);
    ip_id1_id2 = gram_matrix(id1,id2);
    
    % Handle the merge of existing samples
    gram_matrix(:,id1) = alpha1*gram_matrix(:,id1) + alpha2*gram_matrix(:,id2);
    gram_matrix(id1,:) = gram_matrix(:,id1);
    gram_matrix(id1, id1) = alpha1^2*norm_id1 + alpha2^2*norm_id2 + 2*alpha1*alpha2*ip_id1_id2;
    
    gram_vector(id1) = alpha1*gram_vector(id1) + alpha2*gram_vector(id2);
    
    % Handle the new sample
    gram_matrix(:,id2) = gram_vector;
    gram_matrix(id2,:) = gram_matrix(:,id2);
    gram_matrix(id2, id2) = new_sample_norm;
    
    % Update the distance matrix
    distance_matrix(:,id1) = max(gram_matrix(id1, id1) + diag(gram_matrix) - 2*gram_matrix(:,id1),0);
    distance_matrix(id1,:) = distance_matrix(:,id1) ;
    distance_matrix(id1,id1) = inf;
    
    distance_matrix(:,id2) = max(gram_matrix(id2, id2) + diag(gram_matrix) - 2*gram_matrix(:,id2),0);
    distance_matrix(id2,:) = distance_matrix(:,id2) ;
    distance_matrix(id2,id2) = inf;
end


