function secant_method = secant_root(f,R1,R2,tol)
% This function applies the secant method to the function f, with initial
% guesses of roots R1, R2. The number of iterations (iter) must be
% specified.

x = zeros(1,2);
err = ones(1,2);
x(1) = R1;
x(2) = R2;

while abs(err(end))>tol
    j = 1;
    x_new = (x(end-1)*f(x(end))-x(end)*f(x(end-1)))/(f(x(end))-f(x(end-1)));
    x = [x, x_new];
    err_new = round(x(end)-x(end-1),16);
    err = [err, err_new];
    j = j+1;
end

x1 = x.';
err1 = err.';
xanderr = [x1,err1];
results = array2table(xanderr);
results.Properties.VariableNames = {'Test Roots','Error'};
iter = length(x1)-2;
iter_str = num2str(iter);
secant_method = abs(x1(end));

disp(results)
iter_str = append('secant_root took ',iter_str,' iterations to find the root.');
disp(iter_str)