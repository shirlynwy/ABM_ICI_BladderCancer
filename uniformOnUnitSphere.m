function out = uniformOnUnitSphere()

z = 2*rand()-1; % random value on [-1,1]
theta = 2*pi()*rand(); % random value on [0,2*pi)
r = sqrt(1-z^2);
out = [r*[cos(theta);sin(theta)];z];
