function val = constrainValue(val, minVal, maxVal)
%constrainValue: check a value and keep it in a given interval [min, max]
% Args:
% - val (double): value to check
% - minVal(double): minimal bound
% - maxVal(double): maximal bound
% Return:
% - val (double): constrained value

assert(minVal <= maxVal, 'constrainValue: minVal > maxVal!');

if isnan(val)
    warning('constrainValue: isnan detected!');
    val = 0.0;
elseif isinf(val)
    warning('constrainValue: isinf detected!');
    val = 0.0;
elseif val > maxVal
    val = maxVal;
elseif val < minVal
    val = minVal;
%else
    %nop;
end

end

