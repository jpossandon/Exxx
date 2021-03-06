 tk=51;
 cfg_eeg             = eeg_etParams_Exxx('sujid',sprintf('s%02d',tk),...
            'expfolder','/Users/jossando/trabajo/Exxx/');
        filename                = sprintf('s%02d',tk);
    cfg_eeg                 = eeg_etParams_Exxx(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname','stimlock'); 
load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata               
p.times = [1000 1000]
[trls.Cl,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==1'},p.times);            
[trls.Cr,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==2'},p.times);
[trls.Rll,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==5'},p.times);            
[trls.Rl,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==6'},p.times);
[trls.Lr,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==9'},p.times);            
[trls.Lrr,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==10'},p.times);

trls.left = [trls.Cl;trls.Lr;trls.Rll];
trls.right = [trls.Cr;trls.Lrr;trls.Rl];
trls.all   = [trls.left;trls.right];

p.trls_stim = {'Cl','Cr','Lrr','Lr','Rll','Rl','left','right','all'}

%%
p. bsl      = [-.4 0];
p.rref     = 'yes'
for e = 1:9%length(p.trls_stim)
    [ERP.(p.trls_stim{e})] = getERPsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,'ICAem','no');
end

cfgs            = [];
cfgs.parameter  = 'avg';
cfgs.operation  = 'subtract';
ERP.diffLR      = ft_math(cfgs,ERP.left.ICAem,ERP.right.ICAem);

load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
cfgp.xlim       = [-.6 .6]
figure
ft_multiplotER(cfgp,ERP.Lr.ICAem,ERP.Lrr.ICAem)
figure,
ft_multiplotER(cfgp,ERP.all.ICAem)


%%
% TF

% Analysis parameters
p.times_tflock              = [1000 1000];
p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
p.bsl                       = [-.6 -.3]; 
p.rref                     = 'yes';
p.keep                      = 'no';
p.collim                    = [0 2];
p.cfgTFR.channel            = 'all';	
p.cfgTFR.keeptrials         = 'yes';	                
p.cfgTFR.method             = 'mtmconvol';
% p.cfgTFR.taper              = 'hanning'
% p.cfgTFR.width              = 5; 
p.cfgTFR.pad                = 3;
p.cfgTFR.output             = 'fourier';
p.cfgTFR.keeptapers             = 'yes';

p.cfgTFR.foi                = 4:2:90;	

p.cfgTFR.t_ftimwin          = 4./p.cfgTFR.foi;
p.cfgTFR.tapsmofrq          = 0.75*p.cfgTFR.foi;
%plottp(p.cfgTFR)
p.cfgTFR.toi                = (-p.times_tflock(1):20:p.times_tflock(2))/1000;	
plottp(p.cfgTFR)

for e = 1:9%length(p.trls_stim)
    [ERP.(p.trls_stim{e})] = getTFRsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,'ICAem','no',p.cfgTFR);
end
    %%
save('fivetapers','ERP')
%%
cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'divide';
ERP.diffLR      = ft_math(cfgs,ERP.left.ICAem,ERP.right.ICAem);

%%
% Fieldtrip fast plotting
 
load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
%       cfgp.baseline       = p.bsl ;
%       cfgp.baselinetype   = 'db';
% cfgp.trials     = 51:70
%       cfgp.ylim           = [10 18];
%       cfgp.xlim           = [.2 .4];
% cfgp.zlim           = [.7 1.2]
%    cfgp.zlim           = [-2 2];
%   data =TFRallt_RCsac.(p.analysis_type{at});
  data = ERP.diffLR;
% data.powspLFRallt_LU.ICAemUvsCa.powspctrm)
%%
figure
       ft_multiplotTFR(cfgp,data)

 %%
   cfgp.ylim           = [30 70];
       cfgp.xlim           = [.06 .15];
    cfgp.zlim           = [.8 1.2];
 cfgp.comment = 'no'
     figure,ft_topoplotTFR(cfgp,data)