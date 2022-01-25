data = dlmread('data.txt', ',');

str = [];

str = sprintf('%s\n', '{');
for i=1:100
    str = [str, sprintf('%s', '{')];
    str = [str, sprintf('%d', data((i-1)*8+1))];
    for j=(i-1)*8+2:i*8
        str = [str, sprintf(', %d', data(j))];
    end
    str = [str, sprintf('%s, \n', '}')];
end
str = [str, sprintf('%s\n', '};')]