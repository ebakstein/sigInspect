function annot = sigInspectClassifyCov(signal,fs, method, thr, winLength, winAggregPerc,showPlot)
% sigInspectCov
% 
% implemetation of automatic algorithm for stationary segmentation of MER
% signal
%
% Arguments
% ---------
% signal            signal to be segmentated (single channel)
% fs                sampling frequency in Hz
% method            'cov' for auto-covariance (Falkenberg, Aboy 2003) - default
%                   'swt' for stationary wavelet transform (Guarnizo 2006)
% thr                F-test like critical value (recommended 1.2 for COV method)
% winLength     	length of initial segments in seconds(recommended 0.25); 
%                   has to be at most 1/4 of the signal length)
% winAggregPerc     percent of each second marked as artifact necessary to
%                   claim the whole second artifact (default: winLength)
% showPlot          display figure with segmentation+criterion?
%                   Default:false
% 
% Output
% ------
% annot          	annotation of individual seconds based on largest stationary component (boolean)


if(nargin<3)
    method='cov';
end
if(nargin<4)
    switch(method)
        case 'cov'
            thr = 1.2;
        case 'swt'
            thr = 10;
    end
end
if(nargin<5)
    winLength = .25;
end
if(nargin<6)
    winAggregPerc = winLength;
end
if(nargin<6)
    showPlot=false;
end

if((length(signal)/fs)<4*winLength)
   error('winLength has to be at most 1/4 of the signal length!');
end

%% signal normalization
signal=(signal-mean(signal))/std(signal);

% confirm signal is a row vector
if(size(signal,2)==1) 
    signal=signal';
end

%% signal segmentation
winSamples = fs*winLength;
data=buffer(signal,winSamples);
N=length(signal);
Nseg=floor(N/winSamples);
Nsec=ceil(N/fs);

%% variance calculation 
% we use variance of covariance instead of sample variance, because
% segments aren't indepedent observations
switch(method)
    case 'cov' % Falkenberg, Aboy
        covariances=zeros(size(data,1),Nseg);
        for n=1:Nseg
            temp=xcov(data(:,n),'biased');
            covariances(:,n)=temp(round(end/2):end); % we need just one-sided covariance
        end
        variances=var(covariances); 
    case 'swt' % Guarnizo, Orozco        
        LEVELS=8;
        
        % extend signal so that it is divisible by 2^levels
        addCnt=(2^LEVELS-mod(length(signal),2^LEVELS)); 
%         signal=[signal zeros(zcnt,1)];
        if(addCnt>0)
            signal=wextend('1D','symw',signal,addCnt,'r'); % symmetrical extension at the end
        end
                    
        swc=swt(signal,LEVELS,'haar');
%         to=zeros(1,Nseg-1);
%         S{1}=zeros(Nseg,segLen);
%        clear vars V;
        variances=nan(LEVELS+1,Nseg); % wavelet decomposition levels in rows
        for j=1:LEVELS+1, % levels of SWT
            for k=1:Nseg, % segments
                if k*winSamples < length(signal),
%                     S{j}(k,:)=swc(j,1+(k-1)*segLen:k*segLen);
                    s = swc(j,1+(k-1)*winSamples:k*winSamples);
%                 else % no extensions!
% %                     S{j}(k,:)=[swc(j,1+(k-1)*segLen:end) zeros(1,segLen*k-length(signal))];
%                     s=[swc(j,1+(k-1)*winSamples:end) zeros(1,winSamples*k-length(signal))];
                end
                % variance
                variances(j,k)=var(s);
%                 if k>1,
%                     va=sort(vars{j}(k-1:k));
%                     V(j,k-1)=(va(2)/va(1));
%                 end
            end
        end
    otherwise
        error('unknown method: %s',method)
end
%% Sequential hypotesis testing
% statistic test analogous to F test for determining whether two random
% samples have equal variance

%we need to find vn and vd variances for every two segements ie.
% for n=1:length(variances)-1
%     vn(n)=max(variances([n,n+1]));
%     vd(n)=min(variances([n,n+1]));
% end
nvar=size(variances,2);
crit=zeros(nvar,nvar); %matrix which denotes similarity between two segments (where F<FC)
for r=1:nvar
    for c=1:nvar
        if(r~=c)
            F=0;
            for l=1:size(variances,1) % sum over decomp. levels in swt (1 in case of 'cov')
                vn=max(variances(l,[r,c]));
                vd=min(variances(l,[r,c]));
                % "F test" 
                F=F+vn/vd;
            end
            crit(r,c)=F;
            crit(c,r)=F;
        end
    end
end
similar=crit<thr; % test whether they differ a lot

%% Longest segment
% in order to find the longest segment we need to find components of
% relationship of the "similar" matrix and then choose the largest
% component of relationship

segments=zeros(1,size(variances,2)); %denotes which segment belongs to which component
act_comp=1; %in which component we are now

while(~isempty(find(segments==0,1)))
    [~,open]=min(segments); % there is a 0 in the components matrix, therefore min will return index of this 0 (segment that is not in any component yet)
    closed=[];
    
    while(~isempty(open))
        %set component
        segments(open(1))=act_comp;
        
        % children = segments with link to given segment
        sim_row=similar(open(1),:);
        children=find(sim_row==1);
        for ch=1:length(children)
              % current child not in OPEN          AND current child not in CLOSED
            if(isempty(find(open==children(ch),1)))&&(isempty(find(closed==children(ch),1)))
                % ADD current child to OPEN
                open=[open children(ch)];
            end
        end
        
        %remove actual segment
        closed=[closed open(1)]; %actual segment to closed
        open=open(2:end); %remove actual segment from open
    end
    act_comp=act_comp+1;
end

%choose largest component of relationship
comp_length=zeros(1,length(segments)); %there can be at most "length(segments)" of different components of relationship
for c=1:length(segments) 
    comp_length(segments(c))=comp_length(segments(c))+1;
end
[~, comp]=max(comp_length); %"comp" has the number of the longest component
component=find(segments==comp); 

% %%  time indices (unused)
% isegment=zeros(1,winSamples*length(component));
% %longest segment composition
% %keyboard;
% for n=1:length(component)
%     isegment(((n-1)*winSamples+1):(n*winSamples))=((component(n)-1)*winSamples+1):(component(n)*winSamples);
% end
% isegment=isegment(isegment<=length(signal)); % ???
% isegment=isegment(1:length(signal));

%% aggregate segments to 1s windows
boolsegment=true(1,N); 
for ii=1:length(component)
    boolsegment(((component(ii)-1)*winSamples+1):min((component(ii)*winSamples),N)) = false;
end
    
% downsample to 1s (annotation)
annot=false(Nsec,1);
for seci=1:Nsec
    segInds=1+(seci-1)*fs:seci*fs;
    annot(seci)=mean(boolsegment(segInds))>=winAggregPerc; % if at least 20% of the segment is present, mark it! (does not make much difference due to .5s windows)
end



%% plotting
if(showPlot)          
    % display debug plot
    figure;
   s1= subplot(311); xt=(1:length(signal))/fs;plot(xt,signal,'color',[0,0.447,0.741]); 
    hold on; plot(xt(boolsegment),signal(boolsegment),'r.'); title(sprintf('signal + segmentation (%s), FC=%.1f',method,thr));        
    xlabel('time [s]');
    
    % criterion for neighbor components
    critNeigh = crit(Nseg+1:Nseg+1:Nseg*Nseg-1);     
    %compCrit(comp)=[]; % remove the diagonal value
   s2= subplot(312); 
    xn = (winLength:winLength:Nsec-winLength);
    bar(xn,critNeigh,'faceColor',[0.929,0.694,0.125]); hold on; line(xlim,thr*[1 1],'color','r');    
    title('criterion: distance between adjacent segments'); ylabel('criterion value')
    
% criterion for selected components
    critEd = crit; critEd(eye(size(critEd))==1) = Inf;
    compCrit=min(critEd(:,component),[],2)';
    %compCrit(comp)=[]; % remove the diagonal value
   s3 = subplot(313); 
    xs = -winLength/2+( winLength:winLength:Nsec);
    bar(xs,compCrit,'faceColor',[0,0.447,0.741]); hold on; line(xlim,thr*[1 1],'color','r');
    ylim([0 thr*5]);
    title('criterion: min distance from the max component'); ylabel('criterion value')
    linkaxes([s2 s3],'x')
    
end



end