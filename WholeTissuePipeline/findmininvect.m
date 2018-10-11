function quant = findmininvect(x,method,def_quant,belowquant)
%Function in CycIF_Segmentation_normfit_opencellsizeover5_medianpercell
%shave off excess of cell
x = real(x);
x = x(x>0);

if length(x) >0
[n,h] = ksdensity(x);

[~,maxima_loc,widths,~] = findpeaks(n);
num_peaks = length(maxima_loc);
if method == 1 % ie 'minima'
    if num_peaks > 1
        n = -n(maxima_loc(num_peaks-1):maxima_loc(num_peaks));
        h = h(maxima_loc(num_peaks-1):maxima_loc(num_peaks));
        [minimum,minimum_loc] = findpeaks(n);
        quant = h(minimum_loc);
        quant_y  = -n(minimum_loc);
    else
        quant = quantile(x,def_quant);
        quant = quant/belowquant;
        quant_y = 1;
    end
elseif method == 2 % ie 'width'
    try    
        quant = h(real(floor(maxima_loc(num_peaks)-(widths(num_peaks)))));
    catch
        quant = quantile(x,def_quant);
        quant = quant/belowquant;
    end
elseif method == 3 % ie 'only quantile no fitting'
    quant = quantile(x,def_quant);
    quant = quant/belowquant;
end

else
    quant = NaN;
end