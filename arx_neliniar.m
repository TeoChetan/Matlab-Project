function out=arx_neliniar(d,m)

d_size=size(d,2);
d_power=zeros(1,d_size);

for k=1:m
    for i=1:m
        vect=[k zeros(1,d_size-1)];
        for j=1:length(vect)
            if(j>1)
                vect(j)=i;
            end
            if(sum(vect)<=m)
                Perm = unique(permutari(vect),'rows');
                d_power=unique([d_power; Perm],'rows');
            end
        end
    end
end

k=1;
for i=1:size(d_power,1)
    phi=1;
    for j=1:size(d_power,2)
        phi=phi.*d(:,j).^(d_power(i,j));
    end
    Phi(:,k)=phi;
    k=k+1;
end
    out=Phi;
end