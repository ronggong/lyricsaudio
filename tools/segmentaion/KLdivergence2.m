%%% Testing KL Div for PSM
% OUPUTS
% - obs: observation probabilities for each template (in [0,1[)
%
% INPUTS
% - S:          input audio feature (un-normalized) in column
% - templates:  matrix of row templates
% - vdim:       number of dimensions
% - mode:       KL-mode: 0= template on the left (default) d/ 1= template on the right / 2=
% symmetrized
% - KL_exp:     the coefficient before the exp (-.5 by default)
%
% REFERENCES
%               [1] A coupled duration-focused architecture for realtime
%               music to score alignment. Arshia Cont, 2010

function obs = KLdivergence2(S, templates,vdim,mode, KL_exp)

if (nargin < 5)
    KL_exp = -.5;   % -.5 is ad hoc
end

if (nargin < 4)
    mode = 0;
end
    
if (nargin < 3)
    vdim = size(S,1);
else
    S = S(1:vdim,:);
end
S = abs(S); % why do two times abs?
Sn = S;
Dqs = zeros(size(templates,1),size(S,2));

%figure, subplot(311),imagesc(log(Sn))
%ylabel(num2str(size(Dqs)))

for i = 1:size(Sn,2)
    if mod(i,10) == 0
        fprintf('%d,%d\n',i,size(templates,1))
    end
    
    if (sum(S(:,i)) <= 0)
        % if a frame energy is zero, its div Dq(:,i) should be 0.
        
        Dqs(:,i)  = 0.;
        disp(['WARNING: observed spectrogram has no energy at timestep ',num2str(i)])
    else
        % normalize spectrum and fill 0 values
        Sn(:,i) = Sn(:,i)./sum(Sn(:,i));
        Sn((Sn(:,i) <= 0.0000002),i)=0.0000001; % prevent null values
        Sn(:,i) = Sn(:,i)./sum(Sn(:,i)); % renormalize
        
        % compute divergences
        for j = 1:size(templates,1)
            
            switch mode
                case 2 % symmetric
                    Dqs(j,i) = .5*sum( (templates(j,:)'-Sn(:,i)) .* log(templates(j,:)'./Sn(:,i)) );
                case 1 % template on the right
                    if ~isequal(size(Sn(:,i)), size(templates(j,:)'))
                        size(Sn(:,i)), size(templates(j,:)')
                        disp 'NOT EQUAL'
                    end
                    Dqs(j,i) = -sum( Sn(:,i) .* log(templates(j,:)'./Sn(:,i)) );
                otherwise
                    % template on the left, the equation (5) in [1] 
               
                    Dqs(j,i) = sum( templates(j,:)' .* log(templates(j,:)'./Sn(:,i)) );
            end
            %Dsym(j,i) = symkl(templates(j,:)',S(:,i));
        end
    end
%     % debug
%     if ~mod(i,200)
%         disp(['observation : time ',num2str(i)])
%     end
end

%subplot(312),imagesc(log(Sn))
%ylabel(num2str(size(Dqs)))


%% Prepare energy correction
%EnergyCorrection_local = ones(size(Sv_mat)) - exp(-1.0*E_local);

% equation (6) in [1]
obs = exp(KL_exp*Dqs);    
    

%subplot(313), imagesc(obs)
%ylabel(num2str(size(obs)))



end


