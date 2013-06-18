%hw2idx2sub Convert index to subscript assignment.

%TODO: array input, array output for better overall efficiency (k computed
%once)
%NOTE: I should had vectorized this, instead of using it in a loop.

function sub = hw2idx2sub(idx, Val)

%Find the subscripts given the Val and index idx
ndx = idx;
k = [1 cumprod(Val(1:end-1))];
sub = zeros(1, length(Val));

for i = length(Val):-1:1
    vi = rem(ndx-1, k(i)) + 1;
    vj = (ndx - vi)/k(i) + 1;
    sub(i) = vj;
    ndx = vi;
end
