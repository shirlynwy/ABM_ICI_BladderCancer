function M = initializeImmune(M)

% if M.setup.NI0 == 0
imm_fn = fieldnames(M.I);
imm_fn = imm_fn(~startsWith(imm_fn,"tumor_"));
for i = 1:numel(imm_fn)
    M.immunes(:,M.I.(imm_fn{i})) = zeros(0,numel(M.I.(imm_fn{i})));
end
% else
%     M = placeImmune(M,M.setup.NI0);
% end

M.NI = size(M.immunes,1);