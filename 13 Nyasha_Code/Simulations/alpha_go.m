function alpha = alpha_go(Pfa_Desired,N)
    syms x;
    summation = 0;
    for k = 0:1:(N/2-1)
        summation = summation + nchoosek((N/2 -1 +k),k)*(2+x)^(-k);
    end

    Pfa = 2*((1+x)^(-N/2) - (2+x)^(-N/2) * summation);
    eqn = Pfa == Pfa_Desired;
    alpha = double(vpasolve(eqn, x, [0 100]));
end