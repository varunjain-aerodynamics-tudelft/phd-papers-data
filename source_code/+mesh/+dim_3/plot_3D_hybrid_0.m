function [] = plot_3D_hybrid_0(rec_ppp,rec_nodes,Kx,Ky,Kz,boundary_elements_ID,pf)

figure
hold on 

% lift right face
for eln2 = 1:Ky*Kz
    sliceid = 1;
    eln = boundary_elements_ID.left(eln2);
    
    plt_rec_ppp = reshape(rec_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec_nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec_nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec_nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_ppp(sliceid,:,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = boundary_elements_ID.rght(eln2);
    
    plt_rec_ppp = reshape(rec_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec_nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec_nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec_nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_yyy(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_zzz(sliceid,:,:),[pf+1,pf+1]),reshape(plt_rec_ppp(sliceid,:,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

% top bottom face 
for eln2 = 1:Kx*Kz
    sliceid = 1;
    eln = boundary_elements_ID.bttm(eln2);
    
    plt_rec_ppp = reshape(rec_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec_nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec_nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec_nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_ppp(:,sliceid,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = boundary_elements_ID.topp(eln2);
    
    plt_rec_ppp = reshape(rec_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec_nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec_nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec_nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_yyy(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_zzz(:,sliceid,:),[pf+1,pf+1]),reshape(plt_rec_ppp(:,sliceid,:),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

% front back face
for eln2 = 1:Kx*Ky
    sliceid = 1;
    eln = boundary_elements_ID.back(eln2);
    
    plt_rec_ppp = reshape(rec_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec_nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec_nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec_nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_ppp(:,:,sliceid),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
    
    sliceid = pf+1;
    eln = boundary_elements_ID.frnt(eln2);
    
    plt_rec_ppp = reshape(rec_ppp(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_xxx = reshape(rec_nodes.x(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_yyy = reshape(rec_nodes.y(:,:,eln), [pf+1 pf+1 pf+1]);
    plt_rec_zzz = reshape(rec_nodes.z(:,:,eln), [pf+1 pf+1 pf+1]);

    surf(reshape(plt_rec_xxx(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_yyy(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_zzz(:,:,sliceid),[pf+1,pf+1]),reshape(plt_rec_ppp(:,:,sliceid),[pf+1,pf+1]));
%     set(h, 'EdgeColor','none')
end

colorbar('TickLabelInterpreter','latex','FontSize',18)
% title('pressure reconstruction')
colormap jet

set(gca,'TickLabelInterpreter','latex','FontSize',18,'Xcolor','k','Ycolor','k','Zcolor','k');

xlim([-0.1 1.1])
ylim([-0.1 1.1])
zlim([-0.1 1.1])

view(140,15)

end

