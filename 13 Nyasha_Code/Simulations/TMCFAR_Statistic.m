%function ind = TMCFAR_Statistic(Pfa_set,N,T1,T2)
    
    T1 = 18;
    T2 = 2;
    N = 24;
    Pfa_set = 1e-3;
    iterations = 500;
    Pfa_set_vector = ones(iterations,1)*Pfa_set;
    alpha_values = 1:iterations;
    Pfa_values = [];

    for T = alpha_values
        sum_part = [];
        for k = 0:T1 
            sum_part = [sum_part;(factorial(T1)*(-1)^(T1-k)/(factorial(k)*factorial(T1-k)))/(((N-k)/(N-T1-T2))+T)];
        end
        sum_part_1 = sum(sum_part);
        %sum_part_1 = 1;
        M_v_1_T = (factorial(N)/(factorial(T1)*factorial(N-T1-1)*(N-T1-T2)))*sum_part_1;
        M = [M_v_1_T];
        for i = 2:N-T1-T2
            a_i = (N-T1-i+1)/(N-T1-T2-i+1);
            M_v_i_T = a_i/(a_i+T);
            M =[M;M_v_i_T];
        end
        Pfa = prod(M);
        Pfa_values = [Pfa_values;Pfa];
    end
    error = abs(Pfa_set_vector-Pfa_values);
    [val,ind] = min(error);

%end
