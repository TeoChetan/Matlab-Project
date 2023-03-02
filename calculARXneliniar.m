function [out1,out2,out3,out4] = calculARXneliniar(id,val,na,nb,nk,m)
N_id=length(id.y);
N_val = length(val.y);
d=zeros(N_id,na+nb);
d_val = zeros(N_val,na+nb);
% predictie identificare+validare
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


% simulare validare

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

%simulare identificare

y_out_sim_id = zeros(N_id,1);
y_out_sim_id(1) = 0.08;
for k = 2:N_val
    d_temp = zeros(1,na+nb);
    for j = 1:na
        if(k-j <=0)
            d_temp(j)=0;
        else
            d_temp(j) = -y_out_sim_id(k-j);
        end
    end
    for j = na:nb+na
        if(k-j+na <=0)
            d_temp(j) = 0;
        else
            d_temp(j) = id.u(k-j+na);
        end
    end
    y_out_sim_id(k) = arx_neliniar(d_temp,m)*theta;
end

out1 = y_hat_val;
out2 = y_out_sim_val;
out3 = y_hat_id;
out4 = y_out_sim_id;
end