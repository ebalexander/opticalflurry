% analyze images from self-printed mask
function analyzeims
try
    % load images
    load('dataMar29'); % contains cornframes and centframes

    % scene params
    transdist = 0:.5:50; % image locations in translation stage mm
    transdistf = 17; % translation stage coordinate (mm) of most-in-focus image, judged by eye
    ds = -175; % mm, somewhat rough estimate
    df = 233; % mm, from thin lens and ds, f = 100 mm
    A = 51; % mm
    prec = 4; % # sigmas in filter
    cE = (ds*(A/prec)/df)^2; % camera constant
    dtrue = df + (transdist - transdistf); % actual distances at which images were taken

    % derivative scales
    timescale = 0:5:50;
    der1scale = 0:5:45;
    der2scale = 0:10:40;
    
   
    % loop over image sets and derivative scales
    for loc = 1:2 % center and corner patches
        if loc == 1
            I = im2double(centframes);
        else
            I = im2double(cornframes);
        end
        
        % precalculate time derivatives
        for dt = 1:length(timescale)
            Dt = onedfirstderiv(timescale(dt));
            It(:,:,:,dt) = imfilter(I,permute(Dt,[3 1 2]),'same');
        end

        for dxx = 1:length(der2scale)
            clear Ixx Iyy
            % make derivative blurs - Szeliski 4.23
            [D2x D2y] = secondderiv(der2scale(dxx));
            % take spatial derivs
            [ho wo zo] = size(I);
            for i = 2:zo-1 % first and last ims never have valid It
                Ixx(:,:,i) = conv2(I(:,:,i),D2x,'same');
                Iyy(:,:,i) = conv2(I(:,:,i),D2y,'same');
            end

            for dx = 1:length(der1scale)
                clear Ix Iy
                % make derivative kernels - Szeliski 4.21
                [D1x,D1y] = firstderiv(der1scale(dx));
                % take spatial derivs
                for i = 2:zo-1 % first and last ims never have valid It
                    Ix(:,:,i) = conv2(I(:,:,i),D1x,'same');
                    Iy(:,:,i) = conv2(I(:,:,i),D1y,'same');
                end
                
                % trim edge pixels for convolutions
                xtrim = max(length(D1x),length(D2x));
                ytrim = max(length(D1y),length(D2y));
                trim = @(I) I(1+xtrim:end-xtrim,1+ytrim:end-ytrim,1:zo-1);
                Ix = trim(Ix); Ixx = trim(Ixx); Iy = trim(Iy); Iyy = trim(Iyy);
                
                % xIx and yIy
%                 [h w z] = size(Ix);
%                 x = -w/2:w/2-1; y = -h/2:h/2-1;
                % camera calibration: move central pixel
%                 cx = 599.5; cy = 959.5; % central pixel from calibration
%                 offx = floor(ho/2-cx); offy = floor(wo/2-cy);
                collcentx = 960/2; % center pixels of collected data
                collcenty = 600/2;
                if loc == 1 % center
                    collx = 330:520; % collected data in dumb region
                    colly = 270:450;
                else % corner
                    collx = 1:100;
                    colly = 1:100;
                end
                x = collx-collcentx; y = colly-collcenty;
                offx = 0; offy = 0; % assume lens tube alignment
                [X Y] = meshgrid(x-offx,y-offy,ones(1,zo));
                X = trim(X); Y = trim(Y);
                xIx = X.*Ix; yIy = Y.*Iy;

                for dt = 1:length(timescale)
%                     dxx, dx, dt
                    % take appropriately-sized time derivative
                    It = trim(squeeze(It(:,:,:,dt)));
                    % solve equation for u (eq 44) over different It
                    % thresholds
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
                    dtlim = max(1,floor(ceil(7*timescale(dt))/2)); % half length of time deriv filter
                    for i = 1+dtlim:zo-dtlim
                        ix = Ix(:,:,i); iy = Iy(:,:,i); xix = xIx(:,:,i-1); yiy = yIy(:,:,i-1); ixx = Ixx(:,:,i-1); iyy = Iyy(:,:,i-1); it = It(:,:,i);
                        spacederivs = [ix(:) iy(:) xix(:)+yiy(:) ixx(:)+iyy(:)]'*[ix(:) iy(:) xix(:)+yiy(:) ixx(:)+iyy(:)];
                        timederivs = [ix(:) iy(:) xix(:)+yiy(:) ixx(:)+iyy(:)]'*[it(:)];
                        u(:,i) = spacederivs\timederivs;
                    end

                    % solve equation for d (eq 41)
                    u3 = -u(4,:); u2 = -u(3,:); % bad notation
                    d(:,dt,dx,dxx,loc) = (cE*df*u2./(cE*u2-u3));
                end
                save('trialMar29','d','dtrue')
            end
        end
        notifyemma(['dxx = ' num2str(dxx) '/' num2str(length(der2scale))]);
    end
    notifyemma('loc done');
catch ME
   a = ['line ' num2str(ME.stack.line) ', fn ' ME.stack.name ': ' ME.message]
  notifyemma(['line ' num2str(ME.stack.line) ', fn ' ME.stack.name ': ' ME.message]); 
end
end

function [D1x,D1y] = firstderiv(sigma) % Szeliski eq 4.21
    if sigma == 0
        D1x = [-.5 0 .5];
        D1y = D1x';
    else
        w = 2*floor(ceil(7*sigma)/2)+1;
        [xx,yy] = meshgrid(-(w-1)/2:(w-1)/2,-(w-1)/2:(w-1)/2);
        D1x = -xx/(sigma^4).*exp(-0.5*(xx.^2+yy.^2)/(sigma^2));
        D1y = -yy/(sigma^4).*exp(-0.5*(xx.^2+yy.^2)/(sigma^2));
        D1x=D1x./sum(D1x(:));
        D1y=D1y./sum(D1y(:));
    end
end

function [D2x,D2y] = secondderiv(sigma) % Szeliksi eq 4.23
    if sigma == 0
        D2x = [.25 0 -.5 0 .25];
        D2y = D2x';
    else
        w = 2*floor(ceil(7*sigma)/2)+1;
        [xx,yy] = meshgrid(-(w-1)/2:(w-1)/2,-(w-1)/2:(w-1)/2);
        D2x = (1-xx.^2/(2*sigma^2))/(sigma^4).*exp(-0.5*(xx.^2+yy.^2)/(sigma^2));
        D2y = (1-yy.^2/(2*sigma^2))/(sigma^4).*exp(-0.5*(xx.^2+yy.^2)/(sigma^2));
        D2x=D2x./sum(D2x(:));
        D2y=D2y./sum(D2y(:));
    end
end

function D1t = onedfirstderiv(sigma)
    if sigma == 0
        D1t = [-.5 0 .5];
    else
         w = 2*floor(ceil(7*sigma)/2)+1;
         xx = -(w-1)/2:(w-1)/2;
         D1t = exp(-xx.^2/sigma^2);
         D1t = D1t/sum(D1t);
    end
end

% % compare to true dists in a meaningful way...
% a = shiftdim(d,1); b = repmat(dtrue(2:end-1)',[1 16 24]);
% percenterr = (a(:)-b(:))./b(:);
% figure; hist(percenterr,100);

% % show results
% figure; hold on;
% % cmap = colormap(jet(length(thresh)));
% % for t = 1:length(thresh)
%     plot(squeeze(d(:,:,:,:,2)))%,'Color',cmap(t,:));
% % end
% plot(dtrue(2:end-1),'ok');
% set(gca,'XTick',1:25:length(dtrue(2:end-1)));
% dists = int2str(dtrue(2:end-1)');
% set(gca,'XTickLabel',dists(2:25:end-1,:));
% xlabel('true distance (mm)');
% ylabel('reconstructed distance');
% % colorbar;
% % 
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
