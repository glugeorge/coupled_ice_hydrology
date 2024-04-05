function b = bed(x,params)
    xsill = x>params.sill_min & x<params.sill_max;
    xdsill = x>=params.sill_max;
    sill_length = params.sill_max-params.sill_min;
    
    b = params.b0 + params.bx.*x;
    
    b(xsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*(x(xsill)-params.sill_min);

    b(xdsill) = params.b0 + (params.bx*params.sill_min) + params.sill_slope.*sill_length + ...
            params.bx*(x(xdsill)-params.sill_max);

end