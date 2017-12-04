function [ terms ] = fetch_vals_at_layer( layer , measure , N )
% layer 1 to N

    t = combnk(1:N,layer);
    
    terms = [];
    
    for k=1:size(t,1)
        index_terms = cumsum(2.^[t(k,:)-1]); % get all indices in their order provided
        index_terms = index_terms(end); % take just the last one for complete tuple
        terms = [ terms measure(index_terms) ];
    end

end