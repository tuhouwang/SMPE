function [zl, cl] = Column(Layers, dep, c)
% Interpolate variable c onto equidistant depth.
% The input dep and c are cells of the same size, while the output 
% zl and cl are long column vectors.
    zl = [];
    for m = 1 : Layers
        if (isempty(zl) == 1)
            zl = dep{m};
            cl = c{m};
        else
            zl = [zl; dep{m}(2:end)];
            cl = [cl;   c{m}(2:end)];
        end
    end
end
