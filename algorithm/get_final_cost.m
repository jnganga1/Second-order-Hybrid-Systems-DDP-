function [lf,lfx,lfxx] = get_final_cost(x,mu,ind)
% ind indicates which mode (in one gait cycle) is active
lf_mat = Phasefinal_cost(x, mu);
lf = lf_mat(ind);

if nargout > 1
    [lfx_mat, lfxx_mat] = lfPhaseInfo(x, mu); % lfx is gradient not jacobian
     lfx = lfx_mat(:,ind);
     lfxx = lfxx_mat(:,:,ind);
end
end

