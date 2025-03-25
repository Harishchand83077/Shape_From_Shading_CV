
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



% read the image 
E = imread('sphere.jpg'); 
% making sure that it is a grayscale image 
if size(E,3) == 3
    E = rgb2gray(E);
end 
E = double(E); 
% normalizing the image to have maximum of one 
E = E ./ max(E(:)); 
% first compute the surface albedo and illumination direction ... 
[albedo,I,slant,tilt] = estimate_albedo_illumination (E);  



% initializations ... 
[M,N] = size(E); 
% surface normals 
p = zeros(M,N); 
q = zeros(M,N); 
% the surface 
Z = zeros(M,N); 
% surface derivatives in x and y directions 
Z_x = zeros(M,N); 
Z_y = zeros(M,N); 
% maximum number of iterations 
maxIter = 200; 
% the normalized illumination direction 
ix = cos(tilt) * tan(slant); 
iy = sin(tilt) * tan(slant); 
for k = 1 : maxIter 
    % using the illumination direction and the currently estimate 
    % surface normals, compute the corresponding reflectance map. 
    % refer to (57) ... 
    R =  (cos(slant) + p .* cos(tilt)*sin(slant)+ q .* sin(tilt)*sin(slant))./sqrt(1 + p.^2 + q.^2); 
     
    % at each iteration, make sure that the reflectance map is positive at 
    % each pixel, set negative values to zero. 
    R = max(0,R); 
     
    % compute our function f which is the deviation of the computed 
    % reflectance map from the original image ... 
    f = E - R; 
     
    % compute the derivative of f with respect to our surface Z ... refer 
    % to (62) 
    df_dZ = (p+q).*(ix*p + iy*q + 1)./(sqrt((1 + p.^2 + q.^2).^3)* sqrt(1 + ix^2 + iy^2))-(ix+iy)./(sqrt(1 + p.^2 + q.^2)* sqrt(1 + ix^2 + iy^2)); 
     
     
    % update our surface ... refer to (61) 
    Z = Z - f./(df_dZ + eps); % to avoid dividing by zero 
     
    % compute the surface derivatives with respect to x and y 
    Z_x(2:M,:) = Z(1:M-1,:); 
    Z_y(:,2:N) = Z(:,1:N-1); 
     
    % using the updated surface, compute new surface normals, refer to (58) 
    % and (59) 
    p = Z - Z_x; 
    q = Z - Z_y;    
end 
% smoothing the recovered surface 
Z = medfilt2(abs(Z),[21 21]); 
  
% visualizing the result 
figure; 
surfl(Z); 
shading interp; 
colormap gray(256); 
lighting phong; 
 


function [albedo, I, slant, tilt] = estimate_albedo_illumination(E)
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
