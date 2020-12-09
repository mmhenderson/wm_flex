function dprime = get_dprime(predlabs,reallabs,un)
% calculate d' for predicted and actual values
% works for multiple classes
% un is the list of classes in the actual data set


if length(predlabs)~=length(reallabs)
    error('real and predicted labels do not match')
end

hrz=zeros(length(un),1);
fpz=zeros(length(un),1);

nTrials = length(predlabs);



%loop over class labels, get a hit rate and false pos for each (treating
%any other category as non-hit)
for ii=1:length(un)
    
    if sum(reallabs==un(ii))==0 || sum(reallabs~=un(ii))==0
        
        % if one of the categories is completely absent - this will return a
        % nan dprime value
        hr=0.5;
        fp=0.5;
        dprime=nan;
        
        return
        
    else
    
        hr = sum(predlabs==un(ii) & reallabs==un(ii))/sum(reallabs==un(ii));
        fp = sum(predlabs==un(ii) & reallabs~=un(ii))/sum(reallabs~=un(ii));    

        % make sure this never ends up infinite
        % correction from Macmillan & Creelman, use 1-1/2N or 1/2N in place
        % of 1 or 0 
        if hr==0
            hr=1/(2*nTrials);
        end
        if fp==0
            fp=1/(2*nTrials);
        end
        if hr==1
            hr=1-1/(2*nTrials);
        end
        if fp==1
            fp=1-1/(2*nTrials);
        end
        
    end

    %convert to z score (norm dist centered at .50)
    hrz(ii)=norminv(hr,0,1);
    fpz(ii)=norminv(fp,0,1);

end

% dprime is the mean of individual dprimes (for two classes, they will be
% same value)
dprime = mean(hrz-fpz);

end

