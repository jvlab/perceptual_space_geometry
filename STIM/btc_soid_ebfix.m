function [ed,eb_avail,eb_msg]=btc_soid_ebfix(edin,opts)
% [ed,eb_avail,eb_msg]=btc_soid_ebfix(edin,opts) checks whether error bars are sane, and fixes those that are not
% edin: an edirs structure, with fields thresh_mags, thresh_mags_eblo, thresh_mags_ebhi
%
% opts.eb_max:  maximum absolute value of a non-exceptional error bar 
% opts.eb_fill: how to fill in exceptional error bars, defaults to 1
%     0->do not fill
%     n->fill with a fraction of the threshold that is equal to the mean of
%     the n largest fractions (Inf means use all valid error bars in that plane)
%
% ed: an edirs structure with the error bars fixed
% eb_avail: 1 if error bars are available, 0 if not
% eb_msg: a strvcat of messages about what was done
%
%    See also:  BTC_SOID_GETDATA, BTC_SOID_XLSREAD, BC_SOIDFG_VALIDATE.
%
if (nargin<=1) opts=[]; end
opts=filldefault(opts,'eb_fill',1);
opts=filldefault(opts,'eb_max',1); %maximum size of error bar
eb_avail=1;
eb_msg=[];
planes=fieldnames(edin);
for iplane=1:length(planes)
    sane=1;
    plane=planes{iplane};
    d=edin.(plane);
    ndirs=size(d.thresh_mags,1);
    goodeb=[1:ndirs];
    % exceptional eb's are signaled by 0 or by being equal to threshold or by being greater than eb_max in absolute value
    eb_zero=union(find(d.thresh_mags_eblo==0),find(d.thresh_mags_ebhi==0));
    eb_equal=union(find(d.thresh_mags_eblo==d.thresh_mags),find(d.thresh_mags_ebhi==d.thresh_mags));
    eb_max=union(find(abs(d.thresh_mags_eblo)>opts.eb_max),find(abs(d.thresh_mags_ebhi)>opts.eb_max));
    eb_except=union(union(eb_zero,eb_equal),eb_max); %these are the error bars that are not present
    %
    % check sanity of non-exceptional error bars
    msg_except=sprintf('plane %s: %2.0f ebs are exceptional; ',plane,length(eb_except));
    list_avail=setdiff([1:ndirs],eb_except);
    msg_sane=[];
    if (any(d.thresh_mags(list_avail)<=d.thresh_mags_eblo(list_avail)))
        msg_sane='some eblo vals are above thresh; ';
        sane=0;
        eb_avail=0;
    end
    if (any(d.thresh_mags(list_avail)>=d.thresh_mags_ebhi(list_avail)))
        msg_sane=cat(2,msg_sane,'some ebhi vals are below thresh; ');
        sane=0;
        eb_avail=0;
    end
    if isempty(msg_sane)
        msg_sane=sprintf('%2.0f are sane; ',length(list_avail));
    end
    msg_fill=[];
    if ~isempty(eb_except)
        if sane==0
            msg_fill='cannot fill exceptional ebs since sanity check fails';
            eb_avail=0;
        elseif length(list_avail)==0
            msg_fill='cannot fill exceptional ebs since no others are available';
            eb_avail=0;
        elseif opts.eb_fill==0
            msg_fill='cannot fill exceptional ebs since no  fill requested';
            eb_avail=0;
        else
            navg=min(opts.eb_fill,length(list_avail));
            %handle eblo
            frac_lo=(d.thresh_mags(list_avail)-d.thresh_mags_eblo(list_avail))./d.thresh_mags(list_avail);
            frac_lo=flipud(sort(frac_lo(:)));
            mfrac_lo=mean(frac_lo(1:navg));
            d.thresh_mags_eblo(eb_except)=d.thresh_mags(eb_except)*(1-mfrac_lo);
            %handle ebhi
            frac_hi=(-d.thresh_mags(list_avail)+d.thresh_mags_ebhi(list_avail))./d.thresh_mags(list_avail);
            frac_hi=flipud(sort(frac_hi(:)));
            mfrac_hi=mean(frac_hi(1:navg));
            d.thresh_mags_ebhi(eb_except)=d.thresh_mags(eb_except)*(1+mfrac_hi);
            msg_fill=sprintf('exceptional ebs filled with mean of %2.0f others, fracs: [%5.2f %5.2f]',navg,mfrac_lo,mfrac_hi);
        end
    end
    if isempty(msg_fill)
        msg_fill='no ebs to fill in';
    end
    eb_msg=strvcat(eb_msg,cat(2,msg_except,msg_sane,msg_fill));
    ed.(plane)=d;
end
return