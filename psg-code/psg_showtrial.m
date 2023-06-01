%psg_showtrial: show an example trial from the cond (*_sess*.csv) and stimulus (%.png) files
%
% Compares session parameters found with those in the mat file saved by psg_spokes_setup
% Scaling on screen is set here, and independent of what is used during an experiment
%
% First comparison stimulus (2) is to the right of center, then they proceed counterclockwise
%
% 01Jun23:  add domain mode for btc or mpi faces
%
%   See also:  PSG_DEFOPTS, PSG_SPOKES_SETUP, FACES_MPI_PSG_SETUP.
%
comma=',';
if ~exist('stim_angle') stim_angle=2; end %stimulus angle in deg
if ~exist('fix2ctr_angle') fix2ctr_angle=4; end %angle from fixation point to center of stim
if ~exist('disp_angle') disp_angle=12; end %display angle in deg
if ~exist('bkgd_gray') bkgd_gray=128; end %gray level for background
%
domain_mode=getinp('domain: 1->btc, 2->mpi faces','d',[1 2]);
%
pathname=getinp('path','s',[]);
pathname=strrep(pathname,'/',filesep);
pathname=strrep(pathname,'\',filesep);
doscmd=cat(2,'dir/b ',pathname,'\*sess*.csv');
[status,dirlist_string]=dos(doscmd);
ascii_lf=10;
nfiles_dir=length(strfind(dirlist_string,cat(2,'csv',ascii_lf)));
nfiles=0;
dirlist=cell(1,0);
for ifile=1:nfiles_dir
    [dirlist{end+1},dirlist_string]=strtok(dirlist_string);
end
nfiles_avail=length(dirlist);
session_files=sort(dirlist);
nsess_found=length(session_files);
%
mat_files=cell(1,nsess_found);
for isf=1:length(session_files)
    mat_files{isf}=session_files{isf}(1:findstr('_sess',session_files{isf})-1);
    if (domain_mode==2)
        mat_files{isf}=strrep(mat_files{isf},'_gy','');
        mat_files{isf}=strrep(mat_files{isf},'_fc','');
    end
end
mat_files=unique(mat_files);
%
nsess_expected=0;
for imat=1:length(mat_files)
    load(cat(2,pathname,filesep,mat_files{imat}));
    disp(sprintf(' *.mat file %2.0f: %s; contents:',imat,mat_files{imat}));
    disp(s);
    switch domain_mode
        case 1
            nsess_expected=nsess_expected+size(s.sessions,3);
        case 2
            nsess_expected=nsess_expected+2*size(s.sessions,3); %one set of sessions for gray, one for color
    end       
end
disp(sprintf('found %3.0f sessions, %3.0f expected',nsess_found,nsess_expected));
for isess=1:nsess_found
    disp(sprintf(' session file %2.0f->%s',isess,session_files{isess}));
end
sess_choice=getinp('session choice (0 to quit)','d',[0 nsess_found]);
if (sess_choice>0)
    fid=fopen(cat(2,pathname,filesep,session_files{sess_choice}));
    header=fgetl(fid);
    ncompares_found=length(find(header==comma)); %commas between columns
    a=cell(0);
    while 1
        nextline=fgetl(fid);
        if ~ischar(nextline)
            break,
        end
        a{end+1}=nextline;
    end
    fclose(fid);
    ntrials_found=length(a);
    disp(sprintf('this session has %3.0f trials (%3.0f expected), each with one reference and %2.0f comparison stimuli (%3.0f expected)',...
        ntrials_found,size(s.sessions,1),ncompares_found,size(s.sessions,2)-1));
    imgfiles=cell(ntrials_found,1+ncompares_found);
    for itrial=1:ntrials_found
        separators=[0 find(a{itrial}==comma) 1+length(a{itrial})];
        for ifile=1:1+ncompares_found
            imgfiles{itrial,ifile}=a{itrial}(1+separators(ifile):separators(ifile+1)-1);
        end
    end
    ifok=0;
    stims=cell(1,1+ncompares_found);
    while (ifok==0)
        trial_choice=getinp('trial number (0 to quit)','d',[0 ntrials_found]);
        if (trial_choice==0)
            ifok=1;
        else
            %read images
            imsize_found=zeros(1+ncompares_found,2);
            imsize_expected=repmat(s.nchecks*s.nsubsamp,1,2);
            if (domain_mode==1)
                colordepth=1;
            end
            for istim=1:1+ncompares_found
                imgfile=cat(2,imgfiles{trial_choice,istim},'.',s.opts_psg.stim_filetype);
                stims{istim}=double(imread(cat(2,pathname,filesep,imgfile)));
                imsize_found(istim,:)=[size(stims{istim},1),size(stims{istim},2)];
                if (domain_mode==2)
                    stims{istim}=stims{istim}(:,1:(end-mod(end,2)),:); %make sure it's an even height
                end
                if (istim==1) & domain_mode==2
                    imsize_expected(1)=size(stims{istim},1);
                    colordepth=size(stims{istim},3);
                end
                if (istim==1)
                    disp(sprintf('read reference     image %20s, size [%3.0f %3.0f] expected [%3.0f %3.0f]',...
                        imgfile,imsize_found(istim,:),imsize_expected));
                else
                    disp(sprintf('read comparison %2.0f image %20s, size [%3.0f %3.0f] expected [%3.0f %3.0f]',...
                        istim-1,imgfile,imsize_found(istim,:),imsize_expected));
                end
            end
            if ~all(imsize_found==imsize_expected)
                warning('mismatch of image sizes');
            else %display the trial as a figure
                imsize_half=round(imsize_expected/2);
                imrange=cell(1,2);
                for id=1:2
                    imrange{id}=[1:imsize_expected(id)];
                end
                pxls_per_deg=imsize_expected(1)/stim_angle;
                trialfig_half=round(disp_angle*pxls_per_deg/2); %make sure it's an even number
                trialfig_size=2*trialfig_half;
                trialfig=repmat(bkgd_gray,[trialfig_size trialfig_size colordepth]);
                tstring=sprintf('%s: trial %3.0f',session_files{sess_choice},trial_choice);
                %place stimuli in figure
                %first comparison stimulus (2) is to the right of center, then they proceed counterclockwise
                %this assumes matlab standards for imagesc:
                % first dimension is vertical, and with low values at top
                % second dimension is horizontal, with low values at left
                for istim=1:1+ncompares_found
                    if (istim==1)
                        ctr_relpos=[0 0];
                    else
                        ctr_ang=(istim-2)*2*pi/ncompares_found;
                        ctr_relpos=round(fix2ctr_angle*pxls_per_deg*[sin(-ctr_ang),cos(ctr_ang)]);
                    end
                    ctr_pxl=ctr_relpos+trialfig_half-imsize_half;
                    trialfig(ctr_pxl(1)+imrange{1},ctr_pxl(2)+imrange{2},:)=stims{istim};
                end               
                figure;
                set(gcf,'NumberTitle','off');
                set(gcf,'Position',[50 50 1200 800]);
                set(gcf,'Name',tstring);
                imagesc(uint8(trialfig),[0 255]);
                colormap gray
                axis square
                axis off
                %
                axes('Position',[0.02,0.02,0.01,0.01]); %for text
                text(0,0,tstring,'Interpreter','none');
                axis off;
            end
        end
    end
end %sess_choice

