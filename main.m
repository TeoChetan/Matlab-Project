clear all;
close all;
clear console;
clear workspace;
clc;
load('iddata-10.mat');

%afisare date 
plot(id.u), title('Intrare identificare')
figure
plot(id.y),title('Iesire identificare')
figure
plot(val.u), title('Intrare validare')
figure
plot(val.y),title('Iesire validare')

na=1; nb=3;
nk=1;
m=2;
N_id=length(id.y);
N_val = length(val.y);
d=zeros(N_id,na+nb);
d_val = zeros(N_val,na+nb);
% predictie
for k=1:N_id
    for i=1:na
        if(k-i<=0)
            d(k,i)=0;
            d_val(k,i) = 0;
        end
        if(k-i >0 )
            d(k,i)= -id.y(k-i);
            d_val(k,i) = -val.y(k-i);
        end
    end
    for j=na:na+nb+nk-1
        if(k-j<=0)
            d(k,j)=0;
            d_val(k,j)=0;

        end
        if(k-j >0 )
            d(k,j)=id.u(k-j);
            d_val(k,j)=val.u(k-j);
        end
    end
end


phi_id = arx_neliniar(d,m);
phi_val = arx_neliniar(d_val,m);
theta = phi_id\id.y;
y_hat_id = phi_id*theta;
y_hat_val = phi_val*theta;


plot(y_hat_val,'m'),title('Yhat validare')
hold on
plot(val.y,'b')
figure
plot(y_hat_id,'r'),title('Yhat identiifcare')
hold on
plot(id.y,'k')

% simulare


y_out_sim_val = zeros(N_id,1);
y_out_sim_val(1) = 0.08;
for k = 2:N_val
    d_temp = zeros(1,na+nb);
    for j = 1:na
        if(k-j <=0) 
            d_temp(j)=0;
        else
            d_temp(j) = -y_out_sim_val(k-j);
        end
    end
    for j = na:nb+na
        if(k-j+na <=0) 
            d_temp(j) = 0;
        else
            d_temp(j) = val.u(k-j+na);
        end
    end
 y_out_sim_val(k) = arx_neliniar(d_temp,m)*theta;


end

figure
plot(y_out_sim_val,'g'),title('Yhat simulare')
hold on
plot(y_hat_val,'k')

MSE = (1/N_id)*(sum(y_hat_val-y_out_sim_val).^2);



% calcul MSE
M=3;
Na=3; Nb=3; nk=1;
C = nchoosek(1:Na,2);
C = reshape(C(:,perms(1:2)),[],2);
C = [C;[1 1]; [2 2]; [3 3]];
index=1;
MSE_sim_min_val = 1;
MSE_id_min = 1;
MSE_val_min =1;
MSE_sim_min_id = 1;

for m = 1:M
   
    for l = 1:size(C,1)
        na = C(l,1);
        nb = C(l,2);
        [y_hat_val,y_out_sim_val,y_hat_id,y_out_sim_id]=calculARXneliniar(id,val,na,nb,nk,m);
        Mse_sim_val = 0;
        Mse_val = 0;
        Mse_id = 0;
        Mse_sim_id=0;
        for i=1:N_id
            
            Mse_sim_val = Mse_sim_val + ((y_hat_val(i) - y_out_sim_val(i)).^2);
            Mse_id = Mse_id + ((y_hat_id(i) - id.y(i)).^2);
            Mse_val = Mse_val + ((y_hat_val(i) - val.y(i)).^2);
            Mse_sim_id = Mse_sim_id + ((y_hat_id(i) - y_out_sim_id(i)).^2);
            
        end
        MSE_sim_val(index) = (1/N_id)*Mse_sim_val;
        MSE_id(index) = (1/N_id)*Mse_id;
        MSE_val(index) = (1/N_id)*Mse_val;
        MSE_sim_id(index) = (1/N_id)*Mse_sim_id;
       
        if(MSE_sim_val(index)< MSE_sim_min_val)
            MSE_sim_min_val = MSE_sim_val(index);
            na_sim_val_min = na;
            nb_sim_val_min = nb;
            m_sim_val_min = m;
        end
        if(MSE_id(index)< MSE_id_min)
            MSE_id_min = MSE_id(index);
            na_id_min = na;
            nb_id_min = nb;
            m_id_min = m;
        end
        if(MSE_val(index)< MSE_val_min)
            MSE_val_min = MSE_val(index);
            na_val_min = na;
            nb_val_min = nb;
            m_val_min = m;
        end
        if(MSE_sim_id(index)< MSE_sim_min_id)
            MSE_sim_min_id = MSE_sim_id(index);
            na_sim_id_min = na;
            nb_sim_id_min = nb;
            m_sim_id_min = m;
        end
        index = index+1;

    end
   
end

na = [C(:,1);C(:,1);C(:,1)];
nb = [C(:,2);C(:,2);C(:,2)];
m = [ones(size(C,1),1);2*ones(size(C,1),1);3*ones(size(C,1),1)];
MSE_val = MSE_val';
MSE_id = MSE_id';
MSE_sim_val = MSE_sim_val';
MSE_sim_id=MSE_sim_id';
MSEtable = table(m,na,nb,MSE_val,MSE_id,MSE_sim_val,MSE_sim_id);


[y_hat_val,y_out_sim_val,y_hat_id,y_out_sim_id]=calculARXneliniar(id,val,na_sim_val_min,nb_sim_val_min,nk,m_sim_val_min);

figure
plot(y_hat_id,'r'),title('Yhat identificare ideal')
hold on
plot(id.y,'k')
legend('y_h_a_tid','yid')


figure
plot(y_hat_val,'m'),title('Yhat validare ideal')
hold on
plot(val.y,'b')
 legend('y_h_a_tval','yval')


figure
plot(y_out_sim_val,'g'),title('Yhat simulare validare ideal')
hold on
plot(y_hat_val,'k')
 legend('y_s_i_m','y_h_a_tval')


figure
plot(y_out_sim_id,'g'),title('Yhat simulare identificare ideal')
hold on
plot(y_hat_id,'k')
 legend('y_s_i_m','y_h_a_tid')
















