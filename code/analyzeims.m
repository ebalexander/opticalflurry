% analyze images from self-printed mask

% load images
transdist = 1568:5000:511568; % readout from translation stage, in microns
for i = 1:38
    I(:,:,i) = im2double(rgb2gray(imread(['~/opticalflurry/data/trial3/' int2str(transdist(i)) '.tif'])));
end
for i = 40:93
    I(:,:,i) = im2double(rgb2gray(imread(['~/opticalflurry/data/trial3/' int2str(transdist(i)) '.tif'])));
end
for i = 95:length(transdist)
    I(:,:,i) = im2double(rgb2gray(imread(['~/opticalflurry/data/trial3/' int2str(transdist(i)) '.tif'])));
end
I(:,:,39) = (I(:,:,38)+I(:,:,40))/2; % for left out images, take mean
I(:,:,94) = (I(:,:,93)+I(:,:,95))/2;

% scene params
transdistf = 86568; % translation stage coordinate of most-in-focus image, judged by eye
ds = -145; % mm, somewhat rough estimate
df = 340; % mm, from thin lens and ds, f = 100 mm
A = 50; % mm
cEp = (ds*A/df)^2; % partial cE (cE = cEp/Sigma^2, Sigma not really known)
dtrue = df + (transdist - transdistf)/10000;

timescale = 1:50;
der1scale = 1:5:46;
der2scale = 1:10:41;
for dxx = 1:length(der2scale)
    clear Ixx Iyy
    % make derivative blurs
    [g2x g2y] = meshgrid(-6*der2scale(dxx):6*der2scale(dxx),-6*der2scale(dxx):6*der2scale(dxx));
    derivgauss2 = exp(-(g2x.^2+g2y.^2)/(2*der2scale(dxx)^2));
    derivgauss2 = derivgauss2./sqrt(sum(sum(derivgauss2.^2)));
    % make derivative kernels
    deriv2x = conv2(derivgauss2,[.25, 0, -.5, 0, .25]);
    deriv2y = conv2(derivgauss2,[.25, 0, -.5, 0, .25]');
    % take spatial derivs
    [ho wo zo] = size(I);
    for i = 2:zo-1 % first and last ims never have valid It
        Ixx(:,:,i) = conv2(I(:,:,i),deriv2x,'same');
        Iyy(:,:,i) = conv2(I(:,:,i),deriv2y,'same');
    end
    
    for dx = 2:length(der1scale)
        clear Ix Iy
        % make derivative blurs
        [g1x g1y] = meshgrid(-6*der1scale(dx):6*der1scale(dx),-6*der1scale(dx):6*der1scale(dx));
        derivgauss1 = exp(-(g1x.^2+g1y.^2)/(2*der1scale(dx)^2));
        derivgauss1 = derivgauss1./sqrt(sum(sum(derivgauss1.^2)));
        % make derivative kernels
        deriv1x = conv2(derivgauss1,[-.5, 0, .5]);
        deriv2x = conv2(derivgauss2,[.25, 0, -.5, 0, .25]);
        deriv1y = conv2(derivgauss1,[-.5, 0, .5]');
        deriv2y = conv2(derivgauss2,[.25, 0, -.5, 0, .25]');
        % take spatial derivs
        for i = 2:zo-1 % first and last ims never have valid It
            Ix(:,:,i) = conv2(I(:,:,i),deriv1x,'same');
            Iy(:,:,i) = conv2(I(:,:,i),deriv1y,'same');
        end
        % trim edge pixels for convolutions
        xtrim = max(length(deriv1x),length(deriv2x));
        ytrim = max(length(deriv1y),length(deriv2y));
        trim = @(I) I(1+xtrim:end-xtrim,1+ytrim:end-ytrim,:);
        Ix = trim(Ix); Ixx = trim(Ixx); Iy = trim(Iy); Iyy = trim(Iyy);
        % xIx and yIy
        [h w z] = size(Ix);
        x = -w/2:w/2-1; y = -h/2:h/2-1;
        % camera calibration: move central pixel
        cx = 599.5; cy = 959.5; % central pixel from calibration
        offx = floor(ho/2-cx); offy = floor(wo/2-cy);
        [X Y] = meshgrid(x-offx,y-offy,ones(1,z));
        xIx = X.*Ix; yIy = Y.*Iy;

        for dt = 2:length(timescale)
            dxx, dx, dt
            % take time derivative
            It = (circshift(I,[0 0 timescale(dt)]) - ...
                  circshift(I,[0 0 -timescale(dt)]))/(2*timescale(dt)); % could also smooth over time?
            It = trim(It);
            % solve equation for u (eq 44)
%             thresh = 0:max(max(max(abs(It))))/100:max(max(max(abs(It))))*.3;
%             for t = 1:length(thresh)
%                 for i = 1:z
%                     ix = reshape(Ix(:,:,i),[],1);
%                     iy = reshape(Iy(:,:,i),[],1);
%                     xix = reshape(xIx(:,:,i),[],1);
%                     yiy = reshape(yIy(:,:,i),[],1);
%                     ixx = reshape(Ixx(:,:,i),[],1);
%                     iyy = reshape(Iyy(:,:,i),[],1);
%                     it = reshape(It(:,:,i),[],1);
%                     ind = abs(it)>thresh(t);
%                     spacederivs = [ix(ind) iy(ind) xix(ind)+yiy(ind) ixx(ind)+iyy(ind)]'*[ix(ind) iy(ind) xix(ind)+yiy(ind) ixx(ind)+iyy(ind)];
%                     timederivs = [ix(ind) iy(ind) xix(ind)+yiy(ind) ixx(ind)+iyy(ind)]'*[it(ind)];
%                     u(t,i,:) = spacederivs\timederivs;
%                     % u = -u1x -u1y -u2 -u3
%                 end
%             end
            for i = 1+timescale(dt):zo-timescale(dt)
                ix = Ix(:,:,i); iy = Iy(:,:,i); xix = xIx(:,:,i-1); yiy = yIy(:,:,i-1); ixx = Ixx(:,:,i-1); iyy = Iyy(:,:,i-1); it = It(:,:,i);
                spacederivs = [ix(:) iy(:) xix(:)+yiy(:) ixx(:)+iyy(:)]'*[ix(:) iy(:) xix(:)+yiy(:) ixx(:)+iyy(:)];
                timederivs = [ix(:) iy(:) xix(:)+yiy(:) ixx(:)+iyy(:)]'*[it(:)];
                u(:,i) = spacederivs\timederivs;
            end
            
            % solve equation for d (eq 41)
            u3 = -u(4,:); u2 = -u(3,:); % bad notation
            %distcomp = @(sigsq) sum(sum(((cE/sigsq)*df*u2./(cE*u2-u3)-repmat(dtrue,[length(thresh),1])).^2));
            %Sigsq(dt,dx,dxx) = fmincon(distcomp(sigsq),20,-1,0) % honestly kind of a guess, pick best sigma in least squares distance error sense
            sigsq = 4;
            d(:,dt,dx,dxx) = (cEp/sigsq)*df*u2./((cEp/sigsq)*u2-u3);
        end
        save('trial3','d','dtrue')
    end
end
          
% % compare to true dists in a meaningful way...
% a = shiftdim(d,1); b = repmat(dtrue(2:end-1)',[1 16 24]);
% percenterr = (a(:)-b(:))./b(:);
% figure; hist(percenterr,100);

% % show results
% % dfromdf = -375:125:625;
% dfromdf = -35:5:0;
% % trued = df + dfromdf*.254; % convert 1/100ths of inches to mm
% trued = df+dfromdf;
% figure; hold on;
% cmap = colormap(jet(length(thresh)));
% for t = 1:length(thresh)
%     plot(d(t,:),'Color',cmap(t,:));
% end
% plot(trued,'ok');
% set(gca,'XTickLabel',int2str(dfromdf'));
% xlabel(['distance (mm) from focal plane (' num2str(df) 'mm)']);
% ylabel('reconstructed distance');
% colorbar;
% 
% % show images of derivatives
% for t = 1:length(dfromdf);
%     h = figure; 
%     subplot(3,4,1); imshow(I(:,:,t),[]); xlabel('I');
%     subplot(3,4,2); imshow(Ix(:,:,t),[]); xlabel('Ix');
%     subplot(3,4,3); imshow(xIx(:,:,t),[]); xlabel('xIx');
%     subplot(3,4,4); imshow(Ixx(:,:,t),[]); xlabel('Ixx');
%     
%     subplot(3,4,5); imshow(It(:,:,t),[]); xlabel('It');
%     subplot(3,4,6); imshow(Iy(:,:,t),[]); xlabel('Iy');
%     subplot(3,4,7); imshow(yIy(:,:,t),[]); xlabel('yIy');
%     subplot(3,4,8); imshow(Iyy(:,:,t),[]); xlabel('Iyy');
%        
%     subplot(3,2,5); plot(d(:,t)); xlabel('reconstructed dist vs It threshold'); ylabel('should match title + 600 mm')
%     subplot(3,4,11); imshow(xIx(:,:,t)+yIy(:,:,t),[]); xlabel('xIx+yIy');
%     subplot(3,4,12); imshow(Ixx(:,:,t)+Iyy(:,:,t),[]); xlabel('Iyy');
%     
%     suptitle([int2str(dfromdf(t)) ' mm from focal plane']);
%     print(h,['~/Documents/derivpics/smallscale_normed/' int2str(dfromdf(t+1))],'-dpng');
% end
