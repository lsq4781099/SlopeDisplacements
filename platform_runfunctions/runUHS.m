function[im,lambda]=runUHS(sys,opt,h,val)

im      = logsp(0.0001,10,25)';
Tdof    = opt.IM;
opt.im  = repmat(im,1,length(Tdof));
h.id    = h.id{val};
h.p     = h.p(val,:);
h.param = h.param;
h.value = h.value(val,:);
lambda  = runlogictree1(sys,opt,h,1:length(val));
lambda  = permute(lambda,[1,2,3,5,4]);



