function [A,b,Aeq,beq] = setup_lincon()
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    
    addAttachedFiles(gcp,{'dhat.m','fd_weights_full.m','Giesekus_2M3D_v002.m','semhat.m','zwgll.m'});
end
