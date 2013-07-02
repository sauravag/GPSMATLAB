%%% function to perform newton raphson estimation %%%%

function [Ecc_k] = newton_raphson(e,M_k);
    NRnext=0;
    NR=1;
    while abs(NRnext-NR)>1e-12;
        NR=NRnext;
        f=NR-e*sin(NR)-M_k;
        f1=1-e*cos(NR);
        f2=e*sin(NR);
        NRnext=NR-(f/(f1-(f2*f/2*f1))); 
    end;
    Ecc_k=NRnext;		% Eccentric anomaly ESTIMATE
end