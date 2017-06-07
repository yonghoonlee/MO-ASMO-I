function f = obj(x,p)
    [M,Fn] = Giesekus_2M3D_v002(x,p);
    f = [M;Fn];
end