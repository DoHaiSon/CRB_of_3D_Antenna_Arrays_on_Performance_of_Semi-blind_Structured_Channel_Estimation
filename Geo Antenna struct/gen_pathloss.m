function path_loss = gen_pathloss(fc, d, n, sigma, N, M, Nt)
    path_loss = zeros(N, M, Nt);
    for nt=1:Nt
        for nn=1:N
            for mm=1:M
                path_loss(nn, mm, nt) = 20 * log10(4*pi*3*10^8 / fc) ...    % FSPL
                    + 10 * n * log10(d) ...     % PLE
                    + lognrnd(0 , db2mag(sigma));       % SF
            end
        end
    end
    path_loss = db2mag(-path_loss);
end