close all

n = 1000;
xMin = -2*pi;           % x domain min
xMax = 2*pi;            % x domain max
xDom = linspace(xMin,xMax,n);      
dx = xDom(2)-xDom(1);

sigma = 1;
y = (1/(sigma*sqrt(2*pi)))*exp(-.5*(xDom/sigma).^2);

sigma = .01;

y2 = .0125*(1/(sigma*sqrt(2*pi)))*exp(-.5*(xDom/sigma).^2);


y3 = conv(y,y2, 'same');
y4 = fftshift( ifft(fft(y).*fft(y2)) );

plot(xDom,y, xDom, y2, xDom, y3, xDom, y4)
