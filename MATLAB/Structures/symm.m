%%Mirrors the vector such that it is symmetric
%%Author: Andrew Fletcher
%%Email: afletche@ucsd.edu
%%@param data the data to be mirrored
%%@return symm_data the mirrored data
%
%Note: I would like to add a data structure input s.t. I can input a
%matrix or vector. Currently only supports row vectors
function symm_data = symm(data)

%preallocation
memory = zeros(1,length(data));

for i = 1:length(data)
    memory(i) = data(end-(i-1));
end
symm_data = [data memory];