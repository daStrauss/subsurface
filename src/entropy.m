
for x = 1:100
    a = load(['mats/tMat', num2str(x)])
    
    c = hist(a.scrt(:),256)/2000;
    etrp(x) = sum(c(c~=0).*log(c(c~=0)));
    
end
