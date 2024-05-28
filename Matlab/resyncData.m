function new_data = resyncData(new_time, data)
    n = length(new_time);
    nRows = size(data, 1);
    nCols = size(data, 2);
    new_data = NaN(n, nCols);
    new_data(:, 1) = new_time;

    assert(n >= nRows, 'resyncData: n < nRows!');

    i = 1;
    j = 1;
    while i <= nRows && j <= n
        t_i = data(i, 1);
        t_j = new_time(j);
        dt = t_i - t_j;
        while t_i < t_j
            j = j + 1;
            if j > n
                break;
            end
            t_j = new_time(j);
        end
        if (j > 1) && (j <= n)
            dt1 = abs(t_i - new_time(j - 1));
            dt2 = abs(t_i - new_time(j));
            if dt1 < dt2
                new_data(j - 1, 2:end) = data(i, 2:end);
            else
                new_data(j, 2:end) = data(i, 2:end);
            end
        end
        i = i + 1;
    end
end