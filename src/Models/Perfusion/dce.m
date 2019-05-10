classdef dce < AbstractModel % Name your Model
%dce: Compute a Ktrans map using Dynamic Contrast Enhanced data
%<a href="matlab: figure, imshow CustomExample.png ;">Pulse Sequence Diagram</a>
%
% Assumptions:
% (1)FILL
% (2) 
%
% Inputs:
%   DCEData     Dynamic Contrast Enhanced data (4D)
%   T1Map      transverse relaxation time
%   MaskAIF     Mask manually drawn on subject vasculature. (A tool is missing to assist user in this process)
%   (B1Map)     excitation (B1+) fieldmap. Used to correct flip angles. (optional)
%   (Mask)      Binary mask to accelerate the fitting (OPTIONAL)
%
% Outputs:
%   CONCData    Dynamic Concentration data (4D)
%   AIFData     Arterial Input Function (AIF)
%   Ktrans      [mL/100g min]
%   Ve          [mL/100g]
%   (Vb)        [mL/100g]
%   res         Fitting residual
%
%
% Protocol:
%	DCEData  [Tf1 Tf2...Tfn] frames times [s]
%
% Options:
%   Concentration equation
%     'linear'      
%   Method          Method chosen to use in order to fit the data
%     'patlak'      LS (Linear Least Squares)
%     'tofts'       Tofts model NLS (Non-Linear Least Squares)
%     'etofts'      Extended Tofts model NLS (Non-Linear Least Squares)
%
% Example of command line usage:
%   For more examples: <a href="matlab: qMRusage(CustomExample);">qMRusage(CustomExample)</a>
%
% Author: Benoît Bourassa-Moreau (2019)
%
% References:
%   Please cite the following if you use this module:
%     FILL
%   In addition to citing the package:
%     Cabana J-F, Gu Y, Boudreau M, Levesque IR, Atchia Y, Sled JG, Narayanan S, Arnold DL, Pike GB, Cohen-Adad J, Duval T, Vuong M-T and Stikov N. (2016), Quantitative magnetization transfer imaging made easy with qMTLab: Software for data simulation, analysis, and visualization. Concepts Magn. Reson.. doi: 10.1002/cmr.a.21357

properties (Hidden=true)
    onlineData_url = 'Check with Developper Guide/Agah to provide online OSF/qMRLab dataset.';%'https://osf.io/cmg9z/download?version=3';
end
    properties
        MRIinputs = {'DCEData','T1Map','MaskAIF','B1Map','Mask'}; % used in the data panel 
        
        % fitting options
        xnames = {'KTrans','Ve','Vb', 'T1'}; % name of the fitted parameters
        voxelwise = 1; % 1--> input data in method 'fit' is 1D (vector). 0--> input data in method 'fit' is 4D.
        st           = [ 0.05	0.2 0.02 1000]; % starting point
        lb            = [  0      0  0 0]; % lower bound
        ub           = [ 100 10 10 5000 ]; % upper bound
        fx            = [ 0       0       0       1 ]; % fix parameters
        
        % Protocol
%         qibaV4_M0 = 1e3*ones(1,qibaV4_nFrames);
        Prot = struct('SeqParam', ...
            struct('Format',{{'FlipAngle(deg)' 'TR(s)'}},...
            'Mat',[25; 0.005]'),...
            'DCEData', ...
            struct('Format',{{'Time(s)'}},...
            'Mat',(0:0.5:330)'),...
            'Constants', ...
            struct('Format',{{'Hematocrit' 'BrainDensity(g/mL)' 'r1(1/smmol)'}},...
            'Mat',[0.45; 1; 4.5]')); %default protocol for QIBA_v4 phantom
        
        aif = [];
        conc = [];
        
        % Model options
        buttons = {'concetration',{'Linear','Non-linear'},'method',{'Patlak','Tofts','Extended Tofts'}}; %selection buttons
        options= struct();
        
    end
    
    methods
        function obj = dce()
            dbstop in CustomExample.m at 66
            obj.options = button2opts(obj.buttons); % converts buttons values to option structure
        end
        function obj = UpdateFields(obj)
            obj.fx(4) = 1;
        end
        function obj = PrecomputeData(obj,data)
            if ~isfield(data, 'B1Map'), data.B1Map = []; end
            
            flipAngle = obj.Prot.SeqParam.Mat(1,1); % degrees
            TR = obj.Prot.SeqParam.Mat(1,2); % ms
            
            relaxivity1 = obj.Prot.Constants.Mat(1,2);
            
            [obj.conc, obj.aif] = ...
                concandaif(data.DCEData, flipAngle, TR, data.T1Map, ...
                data.MaskAIF, relaxivity1, data.B1Map);
        end
%         function Smodel = equation(obj, x)
%             % Compute the Signal Model based on parameters x. 
%             % x can be both a structure (FieldNames based on xnames) or a
%             % vector (same order as xnames).
%             x = struct2mat(x,obj.xnames);
%             
%             R10 = 1/x(4);
%             S0 = 1e3;
%             flipAngle = obj.Prot.SeqParam.Mat(1,1); % degrees
%             TR = obj.Prot.SeqParam.Mat(1,2); % ms
%             
%             Tvec = obj.Prot.AIFAndTiming.Mat(:,1);
%             Hct = obj.Prot.Constants.Mat(1,2);
%             r1 = obj.Prot.Constants.Mat(2,2);
%             
%             KTrans = x(1); 
%             kep = x(1)/x(2);% kep = KTrans/Ve
%             Vb = x(3);
%             
%                         
%             % % --- Convolution method adapted for variable time vectors
%             % Evaluate IRF model for all ca time repeated for each ct time values
%             dt = Tvec(end)-Tvec(end-1);
%             dtVar = [(Tvec(2:end) - Tvec(1:(end-1))); dt];
%             if any(dtVar - dt > 0.001)
%                 warning('Convolution fails with variable time sampling. Although I don''t know exactly why... should be investigated.')
%             end
%             [tauM, timeM] = meshgrid(Tvec, Tvec);
%             timeConvDiffM = timeM - tauM;
%             
%             biexpIrf = exp(-kep*timeConvDiffM); % (-)
%             biexpIrf(timeConvDiffM<0) = 0;
%             
%             conc = biexpIrf*(aif.*dtVar)+Vb*aif;
%             
%             conc = KTrans*conc;
%             
%             % COMPUTE DATA SIGNAL
%             R1 = r1*conc + R10;
%             Smodel = S0*((1-exp(-TR.*R1))./(1-exp(-TR.*R1).*cos(flipAngle))).*sin(flipAngle);
% 
%         end
        
        function FitResults = fit(obj,data)
            %  Fit data using model equation.
            %  data is a structure. FieldNames are based on property
            %  MRIinputs. 
            
            time = obj.Prot.DCEData.Mat(:,1);
            
            if ~isfield(data, 'B1Map'), data.B1Map = []; end
            
            flipAngle = obj.Prot.SeqParam.Mat(1,1); % degrees
            TR = obj.Prot.SeqParam.Mat(1,2); % ms
            
            relaxivity1 = obj.Prot.Constants.Mat(1,2);
            
            [currentConc, ~] = ...
                concandaif(data.DCEData, flipAngle, TR, data.T1Map, ...
                [], relaxivity1, data.B1Map);
            
            brainDensity = obj.Prot.Constants.Mat(1,3);
            [FitResults.KTrans, FitResults.Ve, FitResults.Vb, AIF] = ...
                DCE_OnSPGR(currentConc, time, obj.aif, ...
                obj.st,obj.lb,obj.ub,obj.fx,brainDensity);
            FitResults.T1 = data.T1Map;
        end
        
        
%         function plotModel(obj, FitResults, data)
%             %  Plot the Model and Data.
%             if nargin<2, qMRusage(obj,'plotModel'), FitResults=obj.st; end
%             
%             %Get fitted Model signal
%             Smodel = equation(obj, FitResults);
%             
%             %Get the varying acquisition parameter
%             Tvec = obj.Prot.Data4D.Mat(:,1); 
%             [Tvec,Iorder] = sort(Tvec);
%             
%             % Plot Fitted Model
%             plot(Tvec,Smodel(Iorder),'b-')
%             
%             % Plot Data
%             if exist('data','var')
%                 hold on
%                 plot(Tvec,data.Data4D(Iorder),'r+')
%                 hold off
%             end
%             legend({'Model','Data'})
%         end
%         
%         function FitResults = Sim_Single_Voxel_Curve(obj, x, Opt, display)
%             % Compute Smodel
%             Smodel = equation(obj, x);
%             % add rician noise
%             sigma = max(Smodel)/Opt.SNR;
%             data.Data4D = random('rician',Smodel,sigma);
%             % fit the noisy synthetic data
%             FitResults = fit(obj,data);
%             % plot
%             if display
%                 plotModel(obj, FitResults, data);
%             end
%         end

    end
end