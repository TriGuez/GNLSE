function radius = gaussRadius(x,y,type)

if strcmp(type,'1/e')
    idx = find(y > max(y)/(exp(1)));
elseif strcmp(type,'1/e2')
    idx = find(y > max(y)/(exp(2)));
elseif strcmp(type,'FWHM')
    idx = find(y > max(y)/2);
else
    fprintf(0,'Unknown type\n');
end
a = min(idx);
b = max(idx);

radius = (x(b)-x(a))/2;