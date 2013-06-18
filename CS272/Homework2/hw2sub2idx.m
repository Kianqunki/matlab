%hw2sub2idx Convert subscript assignment to index.

%TODO: array input, array output for better overall efficiency (k computed
%once)
%NOTE: I should had vectorized this instead of using it in a loop

function idx = hw2sub2idx(sub, Val);

%Compute linear index given the Val and assignment sub
ndx = 1;
k = [1 cumprod(Val(1:end-1))];
for i = 1:length(Val),
    ndx = ndx + (sub(i)-1)*k(i);
end

idx = ndx;