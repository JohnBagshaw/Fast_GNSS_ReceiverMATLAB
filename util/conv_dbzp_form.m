%% This function arrange linear data into DBZP form
% where data_in is an 1xN array and output is 2k x (N/k)
% and first half of each contains rearranged data while
% second half contains zeros.

function [ data_out ] = conv_dbzp_form( data_i, ... % input data
                                        k       ... % block size to convert into
                                        )

    data_out = reshape(data_i, k, length(data_i)/k);
    data_out = [data_out; zeros(size(data_out))];

end

