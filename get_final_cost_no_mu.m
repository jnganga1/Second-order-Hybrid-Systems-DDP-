function [lf,lfx,lfxx] = get_final_cost_no_mu(x,ind)
% ind indicates which mode (in one gait cycle) is active
lf_mat = final_cost_no_mu(x);
lf = lf_mat(ind);

if nargout > 1
    [lfx_mat, lfxx_mat] = lfInfo_no_mu(x); % lfx is gradient not jacobian
     lfx = lfx_mat(:,ind);
     lfxx = lfxx_mat(:,:,ind);
end
end

