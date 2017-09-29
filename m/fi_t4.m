function [ val ] = fi_t4( measure , layer , N )

    if( layer == N )
        val = 0;
        return;
    end
    
    v = fi_t3( measure , layer , N );

    Ts = fetch_vals_at_layer( layer , measure , N );
    
    val = 0;
    for i=1:length(Ts)
        val = val + (Ts(i)-v).^2;
    end

    val = (val) / (length(Ts)-1);
    
end