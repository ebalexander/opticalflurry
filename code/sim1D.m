% u = sim1D(d1,vperp,vpar,dnt,A,p,f,fov,npixels)
% simulates defocus and optical flow combination on
% 1-D fronto-parallel planar texture to find depth.
% d1 and d2 are distance from camera to object at times 1 and 2,
% offset the horizontal difference in camera location,
% A the aperture, f the focal length, and imlen the length of the sensor,
% which is assumed to have arbitrarily small pixels
function [u resid] = sim1D(d1,alphadot,offset,dt,A,p,f,ds,fov,npixels,f0,...
                          noisesigma,firstderivsigma,secondderivsigma,trueu)
  
    % Parameters from input
    theta = fov/2;
    imlen = ds*tand(theta);
    d2 = d1+alphadot*dt;
    d0 = d1-alphadot*dt;
    
    % Preliminary calculations
    alpha0 = -d0/ds;
    alpha1 = -d1/ds;
    alpha2 = -d2/ds;
    [g0 r0] = gauss(d0,f,ds,A,p);
    [g1 r1] = gauss(d1,f,ds,A,p);
    [g2 r2] = gauss(d2,f,ds,A,p);
    sigma0 = r0/(abs(alpha0)*p);    
    sigma1 = r1/(abs(alpha1)*p);
    sigma2 = r2/(abs(alpha2)*p);
    B0 = exp(-2*(pi*f0*sigma0*alpha0)^2);
    B1 = exp(-2*(pi*f0*sigma1*alpha1)^2);
    B2 = exp(-2*(pi*f0*sigma2*alpha2)^2);

    % Generate images    
    X = (-imlen+2*imlen/npixels:2*imlen/npixels:imlen)';
    I0_clean = B0*sin(2*pi*f0*(alpha0*X-offset));
    I1_clean = B1*sin(2*pi*f0*(alpha1*X));
    I2_clean = B2*sin(2*pi*f0*(alpha2*X+offset));
    I0 = I0_clean+normrnd(0,noisesigma,size(X)); 
    I1 = I1_clean+normrnd(0,noisesigma,size(X));
    I2 = I2_clean+normrnd(0,noisesigma,size(X));

%     figure; hold on; plot(I0); plot(I1,'g'); plot(I2,'r'); legend('I0','I1','I2');
      
     % Calculate image derivatives
     % firstderivsigma = noisesigma; % may change this later
     if firstderivsigma
        derivgauss1 = normr(exp(-(-6*firstderivsigma:6*firstderivsigma).^2/(2*firstderivsigma^2)));
     else
         derivgauss1 = 1;
     end
     if secondderivsigma
        derivgauss2 = normr(exp(-(-6*secondderivsigma:6*secondderivsigma).^2/(2*secondderivsigma^2)));
     else
         derivgauss2 = 1;
     end
     deriv1 = conv(derivgauss1,[-.5, 0, .5]);
     deriv2 = conv(derivgauss2,[.25, 0, -.5, 0, .25]);
     Ix = conv(I1,deriv1,'same');
     xIx = X.*Ix;
     Ixx = conv(I1,deriv2,'same');
     It = (I2-I0)/(2*dt); % need whole movie to do smooth derivative
     
     % trim edges
     trim = @(v) v(1+max(length(deriv1),length(deriv2)):end-max(length(deriv1),length(deriv2)));
     trim(Ix); trim(xIx); trim(Ixx); trim(It); 
    
    % Solve for geometry
    if A > 0
        halfspacederivs = [.5*Ix'*Ix    xIx'*Ix      Ixx'*Ix; ...
                            0           .5*xIx'*xIx  Ixx'*xIx;...
                            0           0            .5*Ixx'*Ixx];
        spacederivs = halfspacederivs+halfspacederivs';
        timederivs = [Ix'*It; xIx'*It; Ixx'*It];
        u = spacederivs\timederivs; 
        
        % calculate residual
        resid = [Ix xIx Ixx]*trueu-It;
    
%         % minimize
%         umin = fminunc(@(x)sum((x(1)*Ix+x(2)*IxX+x(3)*Ixx+It).^2),u)
    else % pinhole case, can only get TTC and direction b/c no u(3)
        spacederivs = [Ix'*Ix xIx'*Ix; xIx'*Ix; xIx'*xIx];
        timederivs = [Ix'*It; xIx'*It];
        u = spacederivs\timederivs;
%        ttc = 1/u(2);
%        dir = u(1)/u(2);
        resid = [Ix xIx]*trueu(1:2)-It;
    end
    
end


function [g, r] = gauss(do,f,ds,A,p)
    di = f*do/(do-f); % from thin lens btw do and di
    r = abs(di+abs(ds))*A/2*1/abs(di); % circle of confusion in image coordinates
    r = r * do/ds; % circle in world coordinates
    sigma = r/p; % truncation of gaussian kernel
    x = -r:sign(r):r;
    d = mod(r,1)-1;
    if r ~= 0
        g = exp(-(x-d).^2/(2*sigma^2))/sum(exp(-(x-d).^2/(2*sigma^2)));
    else
        g = double(x==d);
    end
%     r = sign(di-ds)*r;
end