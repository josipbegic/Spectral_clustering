function[]=klaster(k)
%n------>broj tocaka
%k------>broj particija
%dimenzija :D


X=imread('mona.jpg');
n=size(X);

HSV=rgb2hsv(X);

F=[n(1)*n(2)]
for i=1:n(1)*n(2)
    F(i)=HSV(i,i,3)+HSV(i,i,3)*HSV(i,i,2)*sin(HSV(i,i,1))+HSV(i,i,3)*HSV(i,i,2)*cos(HSV(i,i,1));
    
end


for i=1:n(1)
    for j=1:n(2)
                          
        W(i,j)=exp(-((F(i)-(F(j)))^2))*exp(-(i*i+j*j))
        
     end
end

sigma=1;

figure(1)
imagesc(X);


%1. korak
D=diag(W*ones(n,1));


%2. korak
%trazimo prvih K eigenvektora od (D^-1/2)W(D^-1/2), nakon toga racunamo Z
[V,S]=eig(sqrt(D)\W/(sqrt (D)));

[d,order] = sort(diag(S),'descend');  %# sort eigenvalues in descending order
V = V(:,order);

Z=sqrt(D)\V(:,1:k);
%provjeri dal se dobije natrag X po formuli s trece stranice 



%3. korak
XZT=diag((diag(Z*Z')).^(-1/2))*Z;

%4. korak
i=randi(n,1);
R(:,1)=XZT(i,:)'    %5x1
c=zeros(n,1);     %100x1
for j=2:k
    c=c+abs(XZT*R(:,j-1));  %100x1
    [Y,I] = min(c);         
    R(:,j)=XZT(I,:)';
end

%5. korak
FIz=0;

p=0;
while(1)
p=p+1 

%6.korak
XT=XZT*R;
XZ=zeros(n,k);
for i=1:n
    for l=1:k
        
        if (XT(i,l)==max(XT(i,:))) 
        XZ(i,l)=1;
        break
        end
    end
end



%7.korak
[U,O,UT]=svd(XZ'*XZT);
FIt=trace(O);
if abs(FIt-FIz)<(1e-2)
    break
end
abs(FIt-FIz);
FIz=FIt;
R=UT*U';

end

abs(FIt-FIz)
figure(3)
imagesc(XZ);
figure (5)



