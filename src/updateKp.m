% Called by minimization functions to update weights
function K_out = updateKp(K_in, B, I, S, p_new, mask_tune)
Kb_eye_pc_new = K_in.Kb_eye_pc;
Kb_eye_pc_new(mask_tune) = p_new;

Kb_eye = Kb_eye_pc_new(1:length(K_in.Kb_eye));
Kb_pc = Kb_eye_pc_new(length(K_in.Kb_eye)+1:end);

K_out = updateK(B, I, S,  Kb_eye, Kb_pc);
end
