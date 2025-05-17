% btc_sum_eccen:  use models saved by btc_soid_demo to summarize rmse, calculate eccentricities, etc.
%
%  calculations are analytic, based on infinitesimal ellipses -- do not take into
%  account the out-of-plane statistics
%
%   See also:  BTC_PRED_DEMO, BTC_PREDICT.
%
if ~exist('dict') dict=btc_define([]); end
nbtc=length(dict.codel); %10
%
if ~exist('results_fn') results_fn='btc_allraysfixedb_mc_10surrs_test'; end
results_fn=getinp('name of results file name from btc_soid_test (e.g., btc_allraysfixedb_xx_yysurrs)','s',[],results_fn);
r=getfield(load(results_fn),'r');
%
cpairs=[1 2; 1 4; 1 6; 1 10;2 3;2 4;2 6; 2 10;4 5; 4 6; 4 10;6 7;6 8; 6 10];
%
nvariants=length(r);
disp('sy: 0 for no symmetrization around origin, 1 to symmetrize');
disp('ax: 1 to uniformize the thresholds on axes, (0 for not, -1 to keep separate within [bc],[de],[tuvw])');
disp('au: 1 to augment coordinates to match stimuli [default], 0 to ignore augmentation)');
for ivariant=1:nvariants
    disp(sprintf(' variant %2.0f->%s',ivariant,r{ivariant}.setup.label))
end
%
if ~exist('ivarlist') ivarlist=[1:nvariants]; end
ivarlist=getinp('variants to run','d',[1 nvariants],ivarlist);
disp(sprintf('                     ......rmse.......'))
disp(sprintf('    model                    median of      mean of data minus fits'))
disp(sprintf('   variant           raw    surrogates      on axis        off axis'))
for ivarptr=1:length(ivarlist)
    ivariant=ivarlist(ivarptr);
    disp(sprintf('%16s    %7.4f    %7.4f      %7.4f         %7.4f',...
        r{ivariant}.setup.label,...
        r{ivariant}.rmse,...
        r{ivariant}.rmse_med_surrogates,...
        r{ivariant}.mean_data_minus_fit_on_axis,...
        r{ivariant}.mean_data_minus_fit_off_axis));
end
for ivarptr=1:length(ivarlist)
    ivariant=ivarlist(ivarptr);
    disp(sprintf('%16s',r{ivariant}.setup.label));
    %
    qf=r{ivariant}.results.qfit;
    disp('           ellipse fit     ..pooling index..');
    disp('  pair    eccentricity     equal thr   45deg');
    for ipair=1:size(cpairs,1)
        inplane=qf(cpairs(ipair,:),cpairs(ipair,:));
        A=inplane(1,1);
        B=inplane(1,2);
        C=inplane(2,2);
        ecc=sqrt(2*sqrt((A-C)^2+4*B^2)/(A+C+sqrt((A-C)^2+4*B^2)));
        %pooling index along line of equal thresholds, approximated by radial-8 design
        %this is the one that best matches poster and paper
        s2=2*sqrt(A*C)/(A+C); %sin(2*theta)
        c2=(C-A)/(A+C);
        poolindex_et=sqrt(((A+C)+(A-C)*c2+2*B*s2)./((A+C)+(A-C)*c2-2*B*s2));
        %pooling index along w=45-deg line
        poolindex_45=sqrt((A+C+2*B)/(A+C-2*B));
        disp(sprintf('%2s          %6.3f         %6.3f     %6.3f',dict.codel(cpairs(ipair,:)),ecc,poolindex_et,poolindex_45));
    end
end
