% optical flurry noise simulation
% Mar 2015/G2 spring

% fixed params
offset = 0;
prec = 6;
fov = 20;
npixels = 500;
dt = 1; % any smaller and it increases noise in It by division -- this is most optimistic
f = 100; ds = -105; df = (ds*f)/(ds-f); % we have a lens whose f is 100mm, set ds to 105mm
f0 = .5;% set for ~10 cycles in in-focus image
ntrials = 1000; % set as outer loop, easy to stop early, save every time
%fstops = [32 22 16 11 8 5.6 4 2.8 2 1.4 1 .7 .5];
A = 50; % this is our actual setup
cE = (abs(ds)/(df*prec)*A).^2; % camera constant

% independent variables
deltad = .01:.01:.1; % seems resonable
%deltad = .001:.001:.01; % too small
%A = abs(f)./fstops; % we can just stick with our A, bigger found to be
%better
percentd = -.05:.005:.05; % percent d = (do-df)/df 
noise = [0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]/255;
derivsig = noise;

% loop
for i = 1:ntrials
    i
    for n = 1:length(noise)
        now = [i, noise(n)]
        for s = 1:length(derivsig)
            for p = 1:length(percentd)
                for dd = 1:length(deltad)
                    alphadot = deltad(dd)/dt; 
                    d1 = (percentd(p)+1)*df;
                    trueu=[offset*ds/d1,alphadot*ds/d1,...
                           alphadot*(ds/d1)*cE*(d1-df)/d1]';
                    [u_deltad_percentd_derivscale_noise_trial(:,dd,p,s,n,i),...
                     resid(dd,p,s,n,i)] = ...
                        sim1D(d1,alphadot,offset,dt,A,prec,f,abs(ds),...
                              fov,npixels,f0,noise(n),derivsig(s),trueu);
                end
            end
        end
    end
    save('../data/noisesimderivsMar4','u_deltad_percentd_derivscale_noise_trial','resid')
end
 save('../data/noisesimderivsMar4','u_deltad_percentd_derivscale_noise_trial','resid','deltad','percentd','derivsig','noise')

% evaluate accuracy
% u1 = squeeze(u_deltad_percentd_aperture_noise_trial(1,:,:,:,:,:));
u2 = squeeze(u_deltad_percentd_aperture_noise_trial_littledd(2,:,:,:,:,:));
u3 = squeeze(u_deltad_percentd_aperture_noise_trial_littledd(3,:,:,:,:,:));
ntrials = size(u3,5);
% cE = (abs(ds)/(df*prec) * A).^2; % for changing A
% cE = repmat(reshape(cE,[1 1 length(cE)]),[length(deltad) length(percentd)
% 1 length(noise) ntrials]); % for changing A
calcdist = (df*cE.*u(2))./(cE.*u(2)-u(3)) % recovered distance (eq 41)
calcalphadot = calcdist*u(2)/ds % recovered perp velocity (eq 42)
calcbetadot = calcdist*u(1)/ds % recovered parallel velocity (eq 43)
truedist = repmat((percentd+1)*df,[length(deltad) 1 length(A) length(noise) ntrials]);
disterr = abs(calcdist-truedist)./truea;
disterrmean = mean(disterr,5);


% break out by delta d
disterrmeantrunc = disterrmean; disterrmeantrunc(disterrmean>1)=1; disterrmeantrunc(isnan(disterrmean))=1;
cmap = colormap(jet(length(deltad)));
figure; hold on;
for dd = 1:length(deltad)
    subplot(2,2,2); hold on; plot(squeeze(mean(mean(disterrmeantrunc(dd,:,:,:),3),4)),'Color',cmap(dd,:));
    xlabel('percent d from focal plane (~50 mm)'); set(gca,'XTick',1:length(percentd)); set(gca,'XTickLabel',num2str(100*percentd'));
end
h = legend(num2str(deltad'));
set(h,'Position',[.1 .65 .2 .2]);
for dd = 1:length(deltad)
    subplot(2,2,3); hold on; plot(squeeze(mean(mean(disterrmeantrunc(dd,:,:,:),2),4)),'Color',cmap(dd,:)); 
    xlabel('aperture (mm)'); set(gca,'XTick',1:length(A)); set(gca,'XTickLabel',num2str((round(A*1000)/1000)'));
    subplot(2,2,4); hold on; plot(squeeze(mean(mean(disterrmeantrunc(dd,:,:,:),2),3)),'Color',cmap(dd,:)); 
    xlabel('noise level'); set(gca,'XTick',1:length(noise)); set(gca,'XTickLabel',num2str((round(noise*1000)/1000)'));
end

% % plot effective noise from sim1Dderivs:
% sig = sqrt(noise); figure; hold on; plot(sig.^2./(4*sqrt(pi)*sig.^3));
%                                     plot(sig.^2./(16*pi*sig.^6),'g');
% xlabel('noise level'); set(gca,'XTick',1:length(noise)); set(gca,'XTickLabel',num2str((round(noise*1000)/1000)'));
% ylabel('effective noise');
% legend('Ix and It','Ixx');
% 
% % find good measurements
% [errsort I] = sort(abs(aerrmean(:)));
% [b c d e] = ind2sub(size(aerrmean),I(1)) % and so on for I(2), I(3), ...
% 
% % look at overall average error for each var
% aerrmeantrunc = aerrmean; aerrmeantrunc(aerrmean>1)=1; aerrmeantrunc(isnan(aerrmean))=1;
% figure; subplot(2,2,1); plot(squeeze(mean(mean(mean(aerrmeantrunc,2),3),4))); xlabel('delta d (mm)');
%         set(gca,'XTick',1:length(deltad));
%         set(gca,'XTickLabel',num2str(deltad'));
%         subplot(2,2,2); plot(squeeze(mean(mean(mean(aerrmeantrunc,1),3),4))); xlabel('percent d from focal plane (~50 mm)');
%         set(gca,'XTick',1:length(percentd));
%         set(gca,'XTickLabel',num2str(100*percentd'));
%         subplot(2,2,3); plot(squeeze(mean(mean(mean(aerrmeantrunc,1),2),4))); xlabel('aperture (mm)');
%         set(gca,'XTick',1:length(A));
%         set(gca,'XTickLabel',num2str((round(A*1000)/1000)'));
%         subplot(2,2,4); plot(squeeze(mean(mean(mean(aerrmeantrunc,1),2),3))); xlabel('noise level');
%         set(gca,'XTick',1:length(noise));
%         set(gca,'XTickLabel',num2str((round(noise*1000)/1000)'));
% 
% % look at hists of bad measurements
% figure; hist(aerrmean(aerrmean<1));
% for lim = .1:.1:.5
%     [Idd Ipd Ia In] = ind2sub(size(aerrmean),find(aerrmean>lim));
%     figure; subplot(2,2,1); hist(Idd); xlabel('delta d (mm)');
%             set(gca,'XTick',1:length(deltad));
%             set(gca,'XTickLabel',num2str(deltad'));
%             subplot(2,2,2); hist(Ipd); xlabel('percent d from focal plane (~50 mm)');
%             set(gca,'XTick',1:length(percentd));
%             set(gca,'XTickLabel',num2str(100*percentd'));
%             subplot(2,2,3); hist(Ia); xlabel('aperture (mm)');
%             set(gca,'XTick',1:length(A));
%             set(gca,'XTickLabel',num2str((round(A*1000)/1000)'));
%             subplot(2,2,4); hist(In); xlabel('noise level');
%             set(gca,'XTick',1:length(noise));
%             set(gca,'XTickLabel',num2str((round(noise*1000)/1000)'));
%             title(num2str(lim));
% end
