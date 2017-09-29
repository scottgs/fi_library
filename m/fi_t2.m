function [ val ] = fi_t2( measure , layer , N )

    Ts1 = fetch_vals_at_layer( layer , measure , N );
    Ts2 = fetch_vals_at_layer( layer-1 , measure , N );
    
    s1 = sum(Ts1) / length(Ts1);
    s2 = sum(Ts2) / length(Ts2);
    
    val = s1 - s2;
    
end