function [Qxx_tilde,Qux_tilde,Quu_tilde] = regularizer_module(Qxx,Qux,Quu,fgxx,fgux,params)
%Module to do regularization at each tie

method = params.regularizationMethod;
x_size = length(Qxx); 
u_size = length(Quu);
switch method 
    case 2
        %Regularize at each time point 
        Qxx_tilde = Qxx + fgxx; Qux_tilde = Qux + fgux; 
        Qxx_tilde = 0.5 *(Qxx_tilde + Qxx_tilde'); Quu_tilde = 0.5 * (Quu + Quu'); 
        %This follows bc of the schur complement characterization
        %of symmetric matrices  
        %https://www.cis.upenn.edu/~jean/schur-comp.pdf
        % for M =[A B; B^T C]
        % if A > 0 then M > 0 iff the schur complement of A in M > 0
        % if C > 0 then M > 0 iff the schur complement of C in M > 0
        % Below, the schur complement of matrix H in Quu below is Vxx 
        regularization = 0; %for each instance
        H = [Qxx_tilde+1e4*eye(x_size) Qux_tilde'; Qux_tilde Quu_tilde];
        r = length(H);
        [~, p] = chol(H-eye(r)*1e-9);
        while p ~= 0 
            regularization = max(regularization*4, 1e-3);
            Qxx_tilde = Qxx_tilde + eye(x_size)*regularization;                    
            Quu_tilde = Quu_tilde + eye(u_size)*regularization;
            H = [Qxx_tilde+1e4*eye(x_size) Qux_tilde'; Qux_tilde Quu_tilde];
            [~, p] = chol(H-eye(r)*1e-9);
        end        
    case 3
        % Regularize at each time point where regularization is the fgxx term
         p=1;  hyperParam = 1; 
         while p ~= 0 
            Qxx_tilde = Qxx + hyperParam * fgxx - log(hyperParam)*eye(x_size);
            Qxx_tilde = 0.5 *(Qxx_tilde + Qxx_tilde');
            Qux_tilde = Qux + hyperParam * fgux; 
            Quu_tilde = Quu + hyperParam * eye(u_size);

            H = [Qxx_tilde Qux_tilde'; Qux_tilde Quu_tilde];
            r = length(H);
            [~, p] = chol(H-eye(r)*1e-9);
            hyperParam = hyperParam * 0.5;
            if hyperParam < 1e-10
                %if it can't find the proper balance. Just regularize each
                %time point 
                regularization = 0; %for each instance
                p=1;
                while p ~= 0 
                    regularization = max(regularization*4, 1e-3);
                    Qxx_tilde = Qxx + fgxx + eye(x_size)*regularization;                    
                    Quu_tilde = Quu + eye(u_size)*regularization;
                    Qux_tilde = Qux + fgux; 
                    H = [Qxx_tilde+1e4*eye(x_size) Qux_tilde'; Qux_tilde Quu_tilde];
                    [~, p] = chol(H-eye(r)*1e-9);
                end                                
%                 warning('HyperParam reduce below precision. That shouldn''t happen') 
            end
         end       
end



end
