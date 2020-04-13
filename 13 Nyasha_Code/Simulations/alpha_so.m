function alpha = alpha_so(Pfa_set,N)
    syms x;
    summation = 0;
    for k = 0:1:(N/2-1)
        summation = summation + nchoosek((N/2 -1 +k),k)*(2+x)^(-k);
    end

    Pfa = 2*(2+x)^(-N/2)*(summation);
    eqn = Pfa == Pfa_set;
    alpha = double(vpasolve(eqn, x, [0 100]));

end