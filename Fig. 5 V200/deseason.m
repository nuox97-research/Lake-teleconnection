function De_Seasonal = deseason(Variable,Time)
% Variable - for each colomn
% deseason
harm_freq = 1:3;
Timestamp = Time(1);

order = 3;
xt = (Time-Timestamp)/365.25;
x_rad = xt * 2 * pi;
nr_indep = order*2 + 2;
Indep = zeros(height(Variable), nr_indep);
Indep(:,1) = ones(height(Variable),1);
Indep(:,2) = x_rad;

i = 3;
for freq = harm_freq
    cos_freq = cos(x_rad * freq);
    sin_freq = sin(x_rad * freq);
    Indep(:,i) = cos_freq;
    i = i + 1;
    Indep(:,i) = sin_freq;
    i = i + 1;
end

for m = 1:size(Variable,2)
    [b,~,~,~,~] = regress(Variable(:,m),Indep);
    fitted = Indep*b;
    Seasonal(:,m) = fitted;
end

De_Seasonal = Variable-Seasonal;