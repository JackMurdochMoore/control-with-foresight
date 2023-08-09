% A function to turn a vector into a string format useful for file names
% 
% Jack Murdoch Moore, June 2022
% 
function str = vector_to_string(vec)
str = num2str(vec(1));
for ii = 2:numel(vec)
    str = [str, ',', num2str(vec(ii))];
end
end