function [graded_nH] = f_graded_multirun()
    
    close;
    clear;
    clc;
    
    N = 8^4;
    addpath('../')
    % % One parameter case
    % % Neo-Hookean model
    % G_range1 = logspace(-6,-1,N);
    % parfor i = 1:N
    %     disp(i/N)
    %     [t,R] = f_imrv2('neohook',1,'g',G_range1(i))
    %     neohook{i} = [t,R];
    % end
    
    % Four parameter case (G1, G2 l1, l2)
    l1_range4 = logspace(-6,-1,N^(1/4));
    l2_range4 = logspace(-6,-1,N^(1/4));
    G_range4 = logspace(-6,-1,N^(1/4));
    G1_range4 = logspace(-6,-1,N^(1/4));
    [l14, l24, G4, G14] = ndgrid(l1_range4,l2_range4,G_range4, ...
        G1_range4);
    gridPoints4 = cell(length(l1_range4), length(l2_range4), length(G_range4), ...
        length(G1_range4));
    
    for i = 1:length(l1_range4)
        for j = 1:length(l2_range4)
            for k = 1:length(G_range4)
                for l = 1:length(G1_range4)
                    gridPoints4{i,j,k,l} = [l14(i,j,k,l), l24(i,j,k,l), ...
                        G4(i,j,k,l), G14(i,j,k,l)];
                end
            end
        end
    end
    cell4param = reshape(gridPoints4,[1,numel(gridPoints4)]);
    parfor i = 1:N
        %percent = i/N
        [t,R] = f_imrv2('graded',1,'l1',cell4param{i}(1),'l2',cell4param{i}(2),...
            'G',cell4param{i}(3)*cell4param{i}(1)/cell4param{i}(2),...
            'G1',cell4param{i}(4)*cell4param{i}(1)/cell4param{i}(2));
        graded_nH{i} = [t,R];
    end
    graded_nH = reshape(graded_nH,[length(l1_range4), length(l2_range4), ...
        length(G_range4), length(G1_range4)]);
    
end
