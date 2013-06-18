function meas = hw1GetMeasure( Phi, ass )
%hw1GetMeasure  Return the measure value of a factor function, given the assignment.
%On input:
%Phi: The cell data structure that contains the value list, number of values of
%variables, and list of variables.
%ass: A vector that contains the given assignment.
%On output:
%meas: The measurement value corresponding to the given assignment of the
%factor function

List=Phi.List;    %List of measure values of the factor function
Val=Phi.Val;     %Number of values for each variable in the factor function

%TODO: this function should be hw2sub2idx

%Compute linear index
ndx = 1;
k = [1 cumprod(Val(1:end-1))];
for i = 1:length(Val),
    ndx = ndx + (ass(i)-1)*k(i);
end

meas = List(ndx);