function [y] = Integrate_NumNum(r_trg,d,yii,varargin)
r_src = [0 0 0]';
if nargin > 3
    r_src = varargin{1};
end
vol = prod(d);
volvol = vol*vol;
    if r_src == r_trg 
        y = yii; 
    else
        y = volvol/(4*pi)/(norm(r_src-r_trg));
    end
end
