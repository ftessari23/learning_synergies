clear, clc, close all

set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',12);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on'); % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%Number of Muscles
N = 10;

%Limb Lenghts
l1 = 0.35; %Upper Arm
l2 = 0.3; %Fore Arm
l3 = 0.2; %Hand

%Initial Arm Configuration
th10 = 30;
th20 = 45;
th30 = 30;

%Direct Kinematics
[xEE,yEE,xW,yW,xE,yE,xS,yS] = dk_3dof_pend(l1,l2,l3,th10,th20,th30);
plot_3d_pend(xEE,yEE,xW,yW,xE,yE,xS,yS);

%Task-space Jacobian
J = jac_3dof_pend(l1,l2,l3,th10,th20,th30);

%Muscle-space Jacobian
[j_musc,rs,ls,re,le,rw,lw] = jac_m_space(l1,l2,l3,th10,th20,th30);


% Forces
Force_Amp = 10; %[N]
min_ang = 0;
max_ang = 2*pi;
res_ang = (2*pi)/16;
n = input('Enter Direction: [1=+x], [2=+y], [3=-x], [4=-y]');
switch n
    case 1
        phi = [0.0*pi];
        str = 'x_pos';
    case 2
        phi = [0.5*pi];
        str = 'y_pos';
    case 3
        phi = [1.0*pi];
        str = 'x_neg';
    case 4
        phi = [1.5*pi];
        str = 'y_neg';
end
Fx = Force_Amp*cos(phi);
Fy = Force_Amp*sin(phi);
F = [Fx;Fy];

% Torques
J_tr = J';
[U,S,V] = svd(J_tr);
tau = J_tr*F; %[Nm]

% Synthesis of Synergies with SGD
Fm0 = zeros(10,1);
err0 = sum((j_musc*Fm0-tau).^2);

exp_dir = 5;

%Number of Replicates
N_rep = 50;

for replicates = 1:N_rep
    %Initial Starting point
    Fm = Fm0;
    minerr = 100;
    epoch = 1;
    Fm_ev = Fm;
    while minerr > 1e-2 %5e-5
        for i = 1:exp_dir                       
            Fm_dir = Fm;
            Fm_in(:,i) = 50*rand(N,1)-25;
            Fm_dir = Fm_dir+Fm_in(:,i);
            
            for j = 1:N
                if Fm_dir(j) <= 0
                    Fm_dir(j) = 0;
                end
            end

            err(epoch,i) = sum((j_musc*Fm_dir-tau).^2);
            clearvars Fm_dir
        end
        exp = rand(1);
        if exp > 0
            [minerr,idx_min] = min(err(epoch,:));
        else
            idx_min = randi(exp_dir);
            minerr = err(epoch,idx_min);
        end

        Fm = Fm+Fm_in(:,idx_min);
        for j = 1:N
            if Fm(j) <= 0
                Fm(j) = 0;
            end
        end

        % Learning Error
        err_ev(epoch) = minerr;
        if (rem(epoch,2000)) == 0
            err_ev(epoch)
        end
        
        %Force Evolution
        Fm_ev(:,epoch) = Fm;

        clearvars idx_min Fm_in 
        epoch = epoch+1;
    end
    err_ev_rep{replicates} = [err0,err_ev];
    Fm_rep(:,replicates) = Fm;

    clearvars err_ev 
    replicates
end

inf = zeros(N,1);
sup = Inf*ones(N,1);
Fm_ls = lsqlin(j_musc,tau,[],[],[],[],inf,sup);

Fm_rep = [Fm_rep,Fm_ls];
% Vecnorm
Fm_rep_norm = vecnorm(Fm_rep);

% DIEM Calculation
maxFF = max(max(Fm_rep));

Fm_rep_01 = Fm_rep./maxFF;

[exp_center,vard,std_one,orth_med,min_DIEM,max_DIEM] = DIEM_Stat(N,1,0,0);

DIEM = getDIEM(Fm_rep_01,Fm_rep_01,1,0,exp_center,vard);

DIEM = DIEM - eye(size(DIEM)).*diag(DIEM);
DIEM(find(DIEM==0)) = NaN;

% Find Min-Median-Max DIEM for generated replicates
minDIEM_data = min(min(DIEM));
[minDIEM_idxR,minDIEM_idxC] = find(DIEM==min(DIEM(:)));

maxDIEM_data = max(max(DIEM));
[maxDIEM_idxR,maxDIEM_idxC] = find(DIEM==max(DIEM(:)));

%Modal Values
vec_DIEM = sort(DIEM(:));
N_DIEM = length(DIEM);
vec_DIEM = vec_DIEM(1:(N_DIEM*N_DIEM-N_DIEM)/2);
[Nf,edges] = histcounts(vec_DIEM);
[maxx,idxmx] = (max(Nf));
modeDIEM_data = (edges(idxmx)+edges(idxmx+1))/2;
clearvars vec_DIEM Nf edges maxx idxmx
[modeDIEM_idxR,modeDIEM_idxC] = find(abs(DIEM-modeDIEM_data)<5e-2,1,'first');

% Plotting
clc

%Add Force Vectores to the Figure
figure(1),
quiver(xEE*ones(size(Fx)),yEE*ones(size(Fy)),Fx,Fy,0.01,'LineWidth',2,'Color','r'), hold on

figure(1),
nexttile(2);
set(gcf,'Color','w')
quiver(0*ones(size(Fx)),0*ones(size(Fy)),Fx,Fy,'LineWidth',2,'Color','r'), hold on
axis equal
xlabel('$F_x$ [N]','FontName','Times New Roman')
ylabel('$F_y$ [N]','FontName','Times New Roman')
grid on
box off

figure(1),
nexttile(4);
set(gcf,'Color','w')
quiver3(0*ones(size(tau(1,:))),0*ones(size(tau(1,:))),0*ones(size(tau(1,:))),tau(1,:),tau(2,:),tau(3,:),'LineWidth',2), hold on
grid on
xlabel('$\tau_S$ [Nm]','FontName','Times New Roman')
ylabel('$\tau_E$ [Nm]','FontName','Times New Roman')
zlabel('$\tau_W$ [Nm]','FontName','Times New Roman')
axis equal
box off

% Plot
figure()
set(gcf,'Color','white'),
for replicate = 1:N_rep
plot((0:1:length(err_ev_rep{replicate})-1),(err_ev_rep{replicate})), hold on
end
xlim([0 50])
xlabel('Iterations')
ylabel('Error')

% Figure 3
xlab = {'$s_f$','$s_e$','$b_f$','$b_e$','$e_f$','$e_e$','$w_f$','$w_e$','$fw_f$','$fw_e$'};

figure()
set(gcf,'Color','white','Units','inches','Position',[1 1 7 4]),
tt = tiledlayout(1,3);
nexttile(),
barh(1:10,[Fm_rep(:,minDIEM_idxC)';Fm_rep(:,minDIEM_idxR)'])
box off
grid off
yticklabels({'s_f','s_e','b_f','b_e','e_f','e_e','w_f','w_e','f_{w_f}','f_{w_e}'}); 
set(gca,'TickLabelInterpreter', 'tex')
set(gca,'FontName','Times New Roman')
title(strcat('$DIEM_{min}$ = ',num2str(round(minDIEM_data,2))),'FontName','Times New Roman')

nexttile(),
barh(1:10,[Fm_rep(:,modeDIEM_idxC)';Fm_rep(:,modeDIEM_idxR)'])
box off
grid off
yticklabels({'s_f','s_e','b_f','b_e','e_f','e_e','w_f','w_e','f_{w_f}','f_{w_e}'}); 
set(gca,'TickLabelInterpreter', 'tex')
set(gca,'FontName','Times New Roman')
title(strcat('$DIEM_{mode}$ = ',num2str(round(modeDIEM_data,2))),'FontName','Times New Roman')

nexttile(),
barh(1:10,[Fm_rep(:,maxDIEM_idxC)';Fm_rep(:,maxDIEM_idxR)'])
box off
grid off
yticklabels({'s_f','s_e','b_f','b_e','e_f','e_e','w_f','w_e','f_{w_f}','f_{w_e}'}); 
set(gca,'TickLabelInterpreter', 'tex')
set(gca,'FontName','Times New Roman')
title(strcat('$DIEM_{max}$ = ',num2str(round(maxDIEM_data,2))),'FontName','Times New Roman')

ylabel(tt,'Muscles','FontName','Times New Roman')
xlabel(tt,'Force [N]','FontName','Times New Roman')



figure(),
set(gcf,'Color','white'),
h = histogram(DIEM(:),'BinWidth',2); hold on
plot(min_DIEM*ones(10,1),linspace(0,max(h.Values),10),'--k','LineWidth',0.5), hold on
plot(orth_med*ones(10,1),linspace(0,max(h.Values),10),'--k','LineWidth',0.5), hold on
plot(zeros(10,1),linspace(0,max(h.Values),10),'--k','LineWidth',0.5), hold on
box off
grid off
ylabel('Frequency','FontName','Times New Roman')
xlabel('DIEM','FontName','Times New Roman')

figure(),
set(gcf,'Color','white'),
bar(Fm_rep_norm), 
box off
grid off
xlabel('Synergies')
ylabel('Muscle Euclidian Norm')
title('$||F_m||_2$')


figure(),
set(gcf,'Color','white'),
for i = 1:(N_rep+1)
    if i < (N_rep+1)
        err_final(i) = err_ev_rep{1,i}(end);
    else
        err_final(i) = sum((j_musc*Fm_rep(:,N_rep+1)-tau).^2);
    end
end
semilogy(err_final), 
box off
grid off
xlabel('Synergies')
ylabel('Final Error')


function [j_musc,rs,ls,re,le,rw,lw] = jac_m_space(l1,l2,l3,th1,th2,th3)
    % 10 Muscles Model
    % 2x shoulder F/E
        rs = 0.05; %[m]
        ls = 0.14; %[m]
        sf = sqrt(ls^2+rs^2-2*ls*rs*sind(th1)); %[m]
        se = sqrt(ls^2+rs^2+2*ls*rs*sind(th1)); %[m]
    % 2x elbow F/E
        re = 0.035; %[m]
        le = 0.1; %[m]
        ef = sqrt(le^2+re^2-2*le*re*sind(th2)); %[m]
        ee = sqrt(le^2+re^2+2*le*re*sind(th2)); %[m]
    % 2x wrist F/E
        rw = 0.025; %[m]
        lw = 0.08; %[m]
        wf = sqrt(lw^2+rw^2-2*lw*rw*sind(th3)); %[m]
        we = sqrt(lw^2+rw^2+2*lw*rw*sind(th3)); %[m]
    % 2x biarticular shoulder-elbow F/E
        bf = sqrt(l1^2+rs^2+re^2-2*rs*re*cosd(th1+th2)-2*l1*rs*sind(th1)-2*l1*re*sind(th2));
        be = sqrt(l1^2+rs^2+re^2-2*rs*re*cosd(th1+th2)+2*l1*rs*sind(th1)+2*l1*re*sind(th2)); 
    % 2x forearm-wrist F/E
        fwf = sqrt(l2^2+rw^2+re^2-2*rw*re*cosd(th3)-2*l2*rw*sind(th3));
        fwe = sqrt(l2^2+rw^2+re^2-2*rw*re*cosd(th3)+2*l2*rw*sind(th3));

        % j_musc = [ls*rs*cosd(th1)/sf,-ls*rs*cosd(th1)/se,(l1*rs*cosd(th1)-rs*re*sind(th1+th2))/bf,(-l1*rs*cosd(th1)-rs*re*sind(th1+th2))/be,0,0,0,0;...
        %           0,0,(l1*re*cosd(th2)-rs*re*sind(th1+th2))/bf,(-l1*re*cosd(th2)-rs*re*sind(th1+th2))/be,le*re*cosd(th2)/ef,-le*re*cosd(th2)/ee,0,0;...
        %           0,0,0,0,0,0,lw*rw*cosd(th3)/wf,-lw*rw*cosd(th3)/we];

        j_musc = [ls*rs*cosd(th1)/sf,-ls*rs*cosd(th1)/se,(l1*rs*cosd(th1)-rs*re*sind(th1+th2))/bf,(-l1*rs*cosd(th1)-rs*re*sind(th1+th2))/be,0,0,0,0,0,0;...
                  0,0,(l1*re*cosd(th2)-rs*re*sind(th1+th2))/bf,(-l1*re*cosd(th2)-rs*re*sind(th1+th2))/be,le*re*cosd(th2)/ef,-le*re*cosd(th2)/ee,0,0,0,0;...
                  0,0,0,0,0,0,lw*rw*cosd(th3)/wf,-lw*rw*cosd(th3)/we,-(rw*re*sind(th3)-l2*rw*cosd(th3))/fwf,-(rw*re*sind(th3)+l2*rw*cosd(th3))/fwe];

        %Plot Muscle Mechanics
        figure(1),
        
        %Plot Attach Point
        phi = 0:pi/30:(2*pi);
        
        plot(rs*cos(phi),rs*sin(phi),'b-','LineWidth',2), hold on
        plot(re*cos(phi)+l1*cosd(th1),re*sin(phi)+l1*sind(th1),'b-','LineWidth',2), hold on
        plot(rw*cos(phi)+l1*cosd(th1)+l2*cosd(th1+th2),rw*sin(phi)+l1*sind(th1)+l2*sind(th1+th2),'b-','LineWidth',2), hold on
        
        %Plot Muscles - Shoulder
        plot([ls*cosd(th1),0],[ls*sind(th1),rs],'r','LineWidth',2), hold on
        plot([ls*cosd(th1),0],[ls*sind(th1),-rs],'r','LineWidth',2), hold on

        %Plot Muscles - Elbow
        plot([(l1-le)*cosd(th1),l1*cosd(th1)+re*cosd(th1+th2+90)],[(l1-le)*sind(th1),l1*sind(th1)+re*sind(th1+th2+90)],'r','LineWidth',2), hold on
        plot([(l1-le)*cosd(th1),l1*cosd(th1)+re*cosd(th1+th2-90)],[(l1-le)*sind(th1),l1*sind(th1)+re*sind(th1+th2-90)],'r','LineWidth',2), hold on

        %Plot Muscles - Elbow
        plot([(l1-le)*cosd(th1),l1*cosd(th1)+re*cosd(th1+th2+90)],[(l1-le)*sind(th1),l1*sind(th1)+re*sind(th1+th2+90)],'r','LineWidth',2), hold on
        plot([(l1-le)*cosd(th1),l1*cosd(th1)+re*cosd(th1+th2-90)],[(l1-le)*sind(th1),l1*sind(th1)+re*sind(th1+th2-90)],'r','LineWidth',2), hold on
        
        %Plot Muscles - Wrist
        plot([(l1)*cosd(th1)+(l2-lw)*cosd(th1+th2),l1*cosd(th1)+l2*cosd(th1+th2)+rw*cosd(th1+th2+th3+90)],[(l1)*sind(th1)+(l2-lw)*sind(th1+th2),l1*sind(th1)+l2*sind(th1+th2)+rw*sind(th1+th2+th3+90)],'r','LineWidth',2), hold on
        plot([(l1)*cosd(th1)+(l2-lw)*cosd(th1+th2),l1*cosd(th1)+l2*cosd(th1+th2)+rw*cosd(th1+th2+th3-90)],[(l1)*sind(th1)+(l2-lw)*sind(th1+th2),l1*sind(th1)+l2*sind(th1+th2)+rw*sind(th1+th2+th3-90)],'r','LineWidth',2), hold on

        %Plot Muscles - Biarticular 
        plot([0,l1*cosd(th1)+re*cosd(th1+th2+90)],[rs,l1*sind(th1)+re*sind(th1+th2+90)],'r','LineWidth',2), hold on
        plot([0,l1*cosd(th1)+re*cosd(th1+th2-90)],[-rs,l1*sind(th1)+re*sind(th1+th2-90)],'r','LineWidth',2), hold on

        %Plot Muscles - Fore-arm Wrist
        plot([l1*cosd(th1)+re*cosd(th1+th2+90),l1*cosd(th1)+l2*cosd(th1+th2)+rw*cosd(th1+th2+th3+90)],[l1*sind(th1)+re*sind(th1+th2+90),l1*sind(th1)+l2*sind(th1+th2)+rw*sind(th1+th2+th3+90)],'r','LineWidth',2), hold on
        plot([l1*cosd(th1)+re*cosd(th1+th2-90),l1*cosd(th1)+l2*cosd(th1+th2)+rw*cosd(th1+th2+th3-90)],[l1*sind(th1)+re*sind(th1+th2-90),l1*sind(th1)+l2*sind(th1+th2)+rw*sind(th1+th2+th3-90)],'r','LineWidth',2), hold on
end


function plot_3d_pend(xEE,yEE,xW,yW,xE,yE,xS,yS)
    figure('Renderer','painters');
    tiledlayout(2,2),

    nexttile(1,[2 1])
    set(gcf,'Color','w');

    plot(xS,yS,'o','MarkerSize',10,'LineWidth',2,'Color','b'), hold on
    plot([xS xE],[yS yE],'-','LineWidth',2,'Color','k'), hold on

    plot(xE,yE,'o','MarkerSize',10,'LineWidth',2,'Color','b'), hold on
    plot([xE xW],[yE yW],'-','LineWidth',2,'Color','k'), hold on

    plot(xW,yW,'o','MarkerSize',10,'LineWidth',2,'Color','b'), hold on
    plot([xW xEE],[yW yEE],'-','LineWidth',2,'Color','k'), hold on

    plot(xEE,yEE,'o','MarkerSize',10,'LineWidth',2,'Color','r'), hold on

    xlabel('x [m]','FontName','Times New Roman');
    ylabel('y [m]','FontName','Times New Roman');
    
    grid on
    box off
    axis equal

end

function [xEE,yEE,xW,yW,xE,yE,xS,yS] = dk_3dof_pend(l1,l2,l3,th1,th2,th3)
    
    xEE = l1*cosd(th1)+l2*cosd(th1+th2)+l3*cosd(th1+th2+th3);
    yEE = l1*sind(th1)+l2*sind(th1+th2)+l3*sind(th1+th2+th3);

    xW = l1*cosd(th1)+l2*cosd(th1+th2);
    yW = l1*sind(th1)+l2*sind(th1+th2);

    xE = l1*cosd(th1);
    yE = l1*sind(th1);

    xS = 0*ones(size(xEE));
    yS = 0*ones(size(xEE));
    
end

function  J = jac_3dof_pend(l1,l2,l3,th1,th2,th3)
    
    J = [-l1*sind(th1)-l2*sind(th1+th2)-l3*sind(th1+th2+th3), -l2*sind(th1+th2)-l3*sind(th1+th2+th3), -l3*sind(th1+th2+th3);...
         +l1*cosd(th1)+l2*cosd(th1+th2)+l3*cosd(th1+th2+th3), +l2*cosd(th1+th2)+l3*cosd(th1+th2+th3), +l3*cosd(th1+th2+th3)];

end