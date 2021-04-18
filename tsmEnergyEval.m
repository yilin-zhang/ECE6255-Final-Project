function [energy_curves, energy_deviations, e] = tsmEnergyEval(x, y, L, s, Fs)
%tsmEnergyEval.m Returns the energy curves and the energy deviations
% Inputs:
%     x: a vector representing the original
%     y: a vector representing the modified signal
%
% Outputs:
%     energy_curves: a cell that contains 6 energy curves:
%                    Tonal Component of input, Tonal Component of output,
%                    Transient Component of input, Transient Component of
%                    output, Noise Component of input, Noise Component of 
%                    output in order.
%     energy_deviations: a cell that contains 3 energy deviation curves:
%                    Tonal Component, Transient Component, and Noise
%                    Component in order.
%     e: a cell that contains 3 mean energy deviations

    X = spectrogram(x, hamming(L), L/4, L);
    X = abs(X) / max(abs(X), [], 'all');
    L_t = round(0.5 * Fs / (L/4));
    L_f = round(200 / Fs * L);
    [EXs, EXt, EXn] = softMask(X, L_t, L_f);
    
    Y = spectrogram(y, hamming(L), L/4, L);
    Y = abs(Y) / max(abs(Y), [], 'all');
    L_t = round(0.5 * Fs / (L/4) * s);
    L_f = round(200 / Fs * L);
    [EYs, EYt, EYn] = softMask(Y, L_t, L_f);
    
    num_interp_points = length(EYs);
    input_scale = 1:length(EXs);
    if num_interp_points <= 1
        interp_scale = 1;
    else
        interp_scale = 1:(length(EXs)-1)/(num_interp_points-1):length(EXs);
    end
    
    EXs = abs(interp1(input_scale, EXs, interp_scale));
    EXt = abs(interp1(input_scale, EXt, interp_scale));
    EXn = abs(interp1(input_scale, EXn, interp_scale));
    
    LXs = 10*log10(EXs);
    LXt = 10*log10(EXt);
    LXn = 10*log10(EXn);
    
    LYs = 10*log10(EYs);
    LYt = 10*log10(EYt);
    LYn = 10*log10(EYn);
    
    energy_deviations_s = LYs - LXs - (mean(LYs) - mean(LXs));
    energy_deviations_t = LYt - LXt - (mean(LYt) - mean(LXt));
    energy_deviations_n = LYn - LXn - (mean(LYn) - mean(LXn));
    
    es = mean(energy_deviations_s .^ 2);
    et = mean(energy_deviations_t .^ 2);
    en = mean(energy_deviations_n .^ 2);
    
    energy_curves = {LXs, LYs, LXt, LYt, LXn, LYn};
    energy_deviations = {energy_deviations_s, energy_deviations_t, energy_deviations_n};
    e = {es, et, en};
end

function [EXs, EXt, EXn] = softMask(X, L_t, L_f)
    Xs = medfilt2(X, [1, L_t]);
    Xt = medfilt2(X, [L_f, 1]);
    Rs = Xs ./ (Xs + Xt);
    Rt = 1 - Rs;
    Rn = 1 - abs(Rs - Rt);
    
    
    Xmask_s = Rs .* X;
    Xmask_t = Rt .* X;
    Xmask_n = Rn .* X;
    
    EXs = sum(Xmask_s .^ 2, 1);
    EXt = sum(Xmask_t .^ 2, 1);
    EXn = sum(Xmask_n .^ 2, 1);
end