
    function [err_eye, err_pc, err_reg] = getErrorVest(I, pX_curr, p1_curr, KF_X, KF_1, type, mask_tune, JL, RL, WL, R)
        % For tuneVestibularSimultaneous
        
        % Get the current frequency response
        [Ehat_freq, Phat_freq, E_gain, E_phase, P_gain, P_phase] = getFreq(KF_X, I, RL.freqs);
        [Ehat_freq1, Phat_freq1, E_gain1, E_phase1, P_gain1, P_phase1] = getFreq(KF_1, I, RL.freqs);
        E_phase1 = wrapTo180(E_phase1-180);
        E_phase = wrapTo180(E_phase-180);
        
        % Frequency error: difference in gain and phase from desired
        if type==1
            err_eye_freq = ((log(E_gain)- log(RL.E_gain_goal1)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(E_phase - RL.E_phase_goal1)).^2);
            err_pc_freq = ((log(P_gain)- log(RL.P_gain_goal1)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(P_phase - RL.P_phase_goal1)).^2);
            
        elseif type==0
            err_eye_freq = ( ( (log(E_gain)-log(E_gain1)) - (log(RL.E_gain_goal0)-log(RL.E_gain_goal1) ) ).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(  (E_phase-E_phase1) - (RL.E_phase_goal0-RL.E_phase_goal1) )).^2 );
            err_pc_freq = ((log(P_gain)- log(RL.P_gain_goal0)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(P_phase - RL.P_phase_goal0)).^2);
            
        elseif type==2
            err_eye_freq = ( ( (log(E_gain)-log(E_gain1)) - (log(RL.E_gain_goal2)-log(RL.E_gain_goal1) ) ).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(  (E_phase-E_phase1) - (RL.E_phase_goal2-RL.E_phase_goal1) )).^2 );
            err_pc_freq = ((log(P_gain)- log(RL.P_gain_goal2)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(P_phase - RL.P_phase_goal2)).^2);
            
        end
        
        
        % Baseline error: JR data in the dark
        SSE = @(xhat, x) cellfun(@(xhat, x) nansum((xhat - x).^2), xhat, x);
        if type == 1
            mask_JR = ~JL.light & (1:length(JL.light))'<=I.nConds_JR;
            [Ehat_JR, Phat_JR] = getModelJR(KF_X, I, JL.head(mask_JR), JL.target(mask_JR), JL.light(mask_JR), JL.sines(mask_JR));
            
            err_eye_JR = SSE(Ehat_JR, JL.E(mask_JR));
            err_pc_JR = SSE(Phat_JR, JL.P(mask_JR));
        else
            % Specific to learning: only look at Ramachandran & Lisberger frequency data
            err_eye_JR = [];
            err_pc_JR = [];
        end
        
        % STEP
        step_ind = strcmp(JL.conds,'step');
        tt_step = JL.tts{step_ind};
        head_step = JL.head{step_ind}/max(JL.head{step_ind});
        mask_step_ss =  tt_step > I.T_step_ss;  % Start comparing 500 ms in
        mask_step_trans =  tt_step > 0 & tt_step < I.T_step_trans; % Transient is the peak prior to 50 ms - increase to
        
        [Ehat_step, Phat_step] = getModelJR(KF_X, I, {head_step}, {zeros(size(tt_step))}, 0, 0);
        err_eye_step = sum( ( -Ehat_step(mask_step_ss) - WL.eye_ss(type+1)).^2 );
        err_pc_step  = sum( ( Phat_step(mask_step_ss) - WL.pc_ss(type+1)).^2 );
        
        % Get the error in the transient
        if type ~= 1
            temp = Phat_step(mask_step_trans);
            pc_trans_curr = temp(find(abs(temp)==max(abs(temp)),1));
            err_pc_step_trans = (pc_trans_curr - WL.pc_trans(type+1)).^2; % Transient = peak firing in 1st 50 ms after step onset
        else
            % Don't penalize the transient if we are before learning
            err_pc_step_trans  = 0;
        end
        weight_freqs = 100000;
        weight_step = 10000; % changes from 1000
        weight_transient = I.fit_transient*10000; % Don't weight this too much - max measurement of Lisberger 1994 II. may exagerate transient in the data
        
        % TODO: check scale
        err_eye = [err_eye_JR;...
            weight_freqs.*err_eye_freq;... % D
            weight_step*err_eye_step];              % Step steady state
        err_pc = [err_pc_JR; ...
            weight_freqs.*err_pc_freq; ...
            weight_step*err_pc_step;...
            weight_transient*err_pc_step_trans]; % TODO check this!
        
        
        % REGULARIZATION ERROR
        % Create combine regularization matrices for both eye and PC params
        Tall0  = blkdiag(R.T0_eye, R.T0_pc);
        Tall1  = blkdiag(R.T1_eye, R.T1_pc);
        Tall2  = blkdiag(R.T2_eye, R.T2_pc);
        
        % Mask out only the params we care about now
        Tall0_mask = Tall0(mask_tune, mask_tune);
        Tall1_mask = Tall1(mask_tune, mask_tune);
        Tall2_mask = Tall2(mask_tune, mask_tune);
        
        if type == 1
            err_reg_all = I.reg_lambda0*(Tall0_mask*p1_curr).^2 +...
                I.reg_lambda1*(Tall1_mask*p1_curr).^2 + I.reg_lambda2*(Tall2_mask*p1_curr).^2;
        else
            % Penalize CHANGE in weights for learning
            err_reg_all = I.reg_lambda0*(Tall0_mask*(pX_curr-p1_curr)).^2 + ...
                I.reg_lambda1*(Tall1_mask*(pX_curr-p1_curr)).^2 + I.reg_lambda2*(Tall2_mask*(pX_curr-p1_curr)).^2;
        end
        err_reg = sum(err_reg_all);
        
    end
