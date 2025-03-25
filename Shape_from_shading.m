% sphere radius ... 
r = 50; 
% lets define the xy domain (pixel grid)... x and y are data grids  
% (matrices) where the our surface Z = Z(x,y) can be calculated at  
% each point (x,y).
% use finner grid to enhance the resolution 
[x,y] = meshgrid(-1.5*r:0.1:1.5*r,-1.5*r:0.1:1.5*r); 
% surface albedo ... 
albedo = 0.5; 
% illumination direction ... 
I = [0.2 0 0.98]'; 
% surface partial derivates at each point in the pixel grid ... 
p = -x./sqrt(r^2-(x.^2 + y.^2)); % in the x direction 
q = -y./sqrt(r^2-(x.^2 + y.^2)); % in the y direction 
% now lets compute the image brightness at each point in the pixel grid ... 
% the reflectance map will be ... 
R = (albedo.*(-I(1).* p - I(2).* q + I(3)))./sqrt(1 + p.^2 + q.^2); 
  
% the points of the background are those who don't satisfy the sphere 
% equation leading to a negative under the square root of the equation of Z 
% = Z(x,y) 
mask = ((r^2 - (x.^2 + y.^2) >= 0)); 
  
% now we can mask out those points 
R = R .* mask; 
  
% converting the reflectance map to image irradiance (by setting negative 
% values to zeros) 
E = max(0,R); 
  
% converting the image irradiance to a gray scale image 
E = E ./ max(E(:)); 
  
% displaying our synthetic sphere image ... 
figure; 
imshow(E); 
axis off; 
  
% saving our image 
imwrite(E,'sphere.jpg'); 

%% 

% read the image 
E = double(imread('sphere.jpg')); 
  
% normalizing the image to have maximum of one 
E = E ./ max(E(:)); 
  
% compute the average of the image brightness 
Mu1 = mean(E(:)); 
  
% compute the average of the image brightness square 
Mu2 = mean(mean(E.^2)); 
% now lets compute the image's spatial gradient in x and y directions 
[Ex,Ey] = gradient(E); 
% normalize the gradients to be unit vectors 
Exy = sqrt(Ex.^2 + Ey.^2); 
nEx = Ex ./(Exy + eps); % to avoid dividing by zero 
nEy = Ey ./(Exy + eps); 
% computing the average of the normalized gradients 
avgEx = mean(Ex(:)); 
avgEy = mean(Ey(:)); 
% now lets estimate the surface albedo 
gamma = sqrt((6 *(pi^2)* Mu2) - (48 * (Mu1^2))); 
albedo = gamma/pi; 
% estimating the slant 
slant = acos((4*Mu1)/gamma); 
% estimating the tilt 
tilt = atan(avgEy/avgEx); 
if tilt < 0  
tilt = tilt + pi;  
end 
% the illumination direction will be ... 
I = [cos(tilt)*sin(slant) sin(tilt)*sin(slant) cos(slant)]; 
% Display computed values
fprintf('Estimated Albedo: %.4f\n', albedo);
fprintf('Estimated Slant: %.4f\n', slant);
fprintf('Estimated Tilt: %.4f\n', tilt);
fprintf('Estimated Illumination Direction: [%.4f %.4f %.4f]\n', I);

% Display and Save Gradient Images

% Display image gradients Ex and Ey
figure;
subplot(1,2,1);
imshow(mat2gray(Ex));
title('Gradient in X Direction');
imwrite(mat2gray(Ex), 'gradient_x.jpg');

subplot(1,2,2);
imshow(mat2gray(Ey));
title('Gradient in Y Direction');
imwrite(mat2gray(Ey), 'gradient_y.jpg');

% Display normalized gradients nEx and nEy
figure;
subplot(1,2,1);
imshow(mat2gray(nEx));
title('Normalized Gradient in X Direction');
imwrite(mat2gray(nEx), 'normalized_gradient_x.jpg');

subplot(1,2,2);
imshow(mat2gray(nEy));
title('Normalized Gradient in Y Direction');
imwrite(mat2gray(nEy), 'normalized_gradient_y.jpg');



E = imread('sphere.jpg'); 
% making sure that it is a grayscale image 
if size(E,3) == 3
    E = rgb2gray(E);
end

% downsampling to speedup 
E = E(1:2:end,1:2:end); 
E = double(E); 
% normalizing the image to have maximum of one 
E = E ./ max(E(:)); 

% first compute the surface albedo and illumination direction ... 
[albedo,I] = estimate_albedo_illumination (E); 

% Initializations ... 
Z = zeros(M,N); 
% assign very bright pixels with large depth value, the negative is used because it 
% was observed that this algorithm always generate concave surfaces (curved  
% inward), this is an ambiguity in shape from shading in general, SFS algorithm can  
% not distinguish between concave and convex surfaces, both are the same. 
Z(find(E>0.75)) = -100 .* E(find(E>0.75)); 
% getting the surface normals from the initial depth map 
[p,q] = gradient(Z); 
[M,N] = size(E); 
% surface normals 
p = zeros(M,N); 
q = zeros(M,N); 
% the second order derivatives of surface normals 
p_ = zeros(M,N); 
q_ = zeros(M,N); 
% the estimated reflectance map 
R = zeros(M,N); 
% the controling parameter 
lamda = 1000; 
% maximum number of iterations 
maxIter = 2000; 
% The filter to be used to get the second order derivates of the surface 
% normals. 
w = 0.25*[0 1 0;1 0 1;0 1 0]; % refer to equations (39) and (40)  
% wx and wy in (47) and (48) 
[x,y] = meshgrid(1:N,1:M); 
wx = (2.* pi .* x) ./ M; 
wy = (2.* pi .* y) ./ N; 

for k = 1 : maxIter 
  % compute the second order derivates of the surface normals 
p_ = conv2(p,w,'same'); 
q_ = conv2(q,w,'same'); 
% Using the computed surface albedo, illumination direction and 
% surface normals, compute estimation for the reflectance map. 
% refer to (16) 
R = (albedo.*(-I(1).* p - I(2).* q + I(3)))./sqrt(1 + p.^2 + q.^2);   
% Compute the partial derivatives of the reflectance map with respect 
% to p and q. it will be the differenation of (16) with respect to p 
% and q 
pq = (1 + p.^2 + q.^2); 
dR_dp =  (-albedo*I(1) ./ (pq .^(1/2))) + (-I(1) * albedo .* p - I(2) * albedo .* q +  I(3) * albedo) .* (-1 .* p .* (pq .^(-3/2))); 

dR_dq =  (-albedo*I(2) ./ (pq .^(1/2))) + (-I(1) * albedo .* p - I(2) * albedo .* q +  I(3) * albedo) .* (-1 .* q .* (pq .^(-3/2))); 
% Compute the newly estimated surface normals ... refer to (41) and 
% (42) 
p = p_ + (1/(4*lamda))*(E-R).*dR_dp; 
q = q_ + (1/(4*lamda))*(E-R).*dR_dq; 
% Compute the Fast Fourier Transform of the surface normals. 
Cp = fft2(p); 
Cq = fft2(q); 
% Compute the Fourier Transform of the surface Z from the Fourier 
% Transform of the surface normals ... refer to (46) ... 
C = -i.*(wx .* Cp + wy .* Cq)./(wx.^2 + wy.^2); 
% Compute the surface Z using inverse Fourier Transform ... refer to 
% (45) 
Z = abs(ifft2(C)); 
% Compute the integrable surface normals .. refer to (47) and (48) 
p = ifft2(i * wx .* C); 
q = ifft2(i * wy .* C); 
% saving intermediate results 
if (mod(k,100) == 0) 
save(['Z_after_' num2str(k) '_iterations.mat'],'Z'); 
end 
end 
% visualizing the result 
figure; 
surfl(Z); 
shading interp; 
colormap gray(256); 
lighting phong; 

function [albedo, I] = estimate_albedo_illumination(E)
    % Compute brightness statistics
    Mu1 = mean(E(:));
    Mu2 = mean(mean(E.^2));

    % Compute estimated albedo
    gamma = sqrt((6 * (pi^2) * Mu2) - (48 * (Mu1^2)));
    albedo = gamma / pi;

    % Estimate slant
    slant = acos((4 * Mu1) / gamma);

    % Compute the image gradients
    [Ex, Ey] = gradient(E);

    % Compute averages of gradients
    avgEx = mean(Ex(:));
    avgEy = mean(Ey(:));

    % Estimate tilt
    tilt = atan(avgEy / avgEx);
    if tilt < 0
        tilt = tilt + pi;
    end

    % Compute illumination direction
    I = [cos(tilt) * sin(slant), sin(tilt) * sin(slant), cos(slant)];
end