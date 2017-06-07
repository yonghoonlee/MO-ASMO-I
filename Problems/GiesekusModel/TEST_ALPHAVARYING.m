testi=3;
testx=DATA{22,9}(testi,:);
j=1;
testf=feval(problem.objfun,testx,problem.p);
testM(j)=testf(1);
testFn(j)=testf(2);
for i = 0.03:0.01:0.07
    j=j+1;
    testx(end-1)=i;
    testf=feval(problem.objfun,testx,problem.p);
    testM(j)=testf(1);
    testFn(j)=testf(2);
end
figure;
plot(DATA{22,11}(:,1),DATA{22,11}(:,2),'bx');
hold on;
plot(DATA{22,11}(testi,1),DATA{22,11}(testi,2),'gs');
plot(testM,testFn,'ro');