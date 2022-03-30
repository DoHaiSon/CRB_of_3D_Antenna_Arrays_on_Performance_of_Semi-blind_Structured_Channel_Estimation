function path_loss = gen_pathloss(alpha, d, beta, sigma, row, col)
    path_loss = zeros(row, col);
    for r=1:row
        for c=1:col
            path_loss(r, c) = 10 * alpha * log10(d) + beta ...
                + normrnd(0, sigma^2);
        end
    end
    path_loss = db2mag(-path_loss);
end

