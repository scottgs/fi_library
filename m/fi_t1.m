function [ val ] = fi_t1( measure , layer , N )

    Ts = fetch_vals_at_layer( layer , measure , N );
    
    val = 1 - ( sum(Ts) / length(Ts) );
    
end