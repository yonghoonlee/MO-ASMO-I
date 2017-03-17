function x = descaling(xs,lb,ub)
    lb = lb';
    ub = ub';
    repeat = size(xs,1);
    lb = repmat(lb,repeat,1);
    ub = repmat(ub,repeat,1);
    x = lb + (ub - lb).*xs;
end