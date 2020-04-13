function alpha = SO_CA_CFAR_Statistic(Pfa_Desired,N)
    syms a;
    
    summation = 0;
    for k = 0:1:(N/2-1)
        summation = summation + nchoosek((N/2 -1 +k),k)*(2+a)^(-k);
    end

    Pfa = 2*(2+a)^(-N/2)*(summation);

    eqn = Pfa == Pfa_Desired;

    alpha = double(vpasolve(eqn, a, [0 1e5]));

end