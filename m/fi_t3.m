function [ val ] = fi_t3( measure , layer , N )

    Ts = fetch_vals_at_layer( layer , measure , N );
    
    val = sum(Ts) / length(Ts);

end