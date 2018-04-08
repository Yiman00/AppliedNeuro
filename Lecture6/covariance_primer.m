% Generate white noise (noise containing many frequencies with equal
% intensities 

n = 500; % n = num of datapoints
a(:,1) = normrnd(0,1,n,1); % mean of 0, std of 1, n x 1 matrix
a(:,2) = normrnd(0,1,n,1); % same parameters to make 2nd dim of matrix a

b(:,1) = normrnd(0,1,n,1); % same parameters to make matrix b
b(:,2) = b(:,1)*0.5 + 0.5*normrnd(0,1,n,1); % correlate 2nd dimension with 1st dimension of matrix b

figure;
subplot(1,2,1); scatter(a(:,1), a(:,2), 'b','.'); title('Uncorrelated Noise'); 
hold on;
subplot(1,2,2); scatter(b(:,1), b(:,2), 'b', '.'); title('Correlated Noise'); 
