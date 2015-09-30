function[]=klaster(n,k)
%n------>broj tocaka
%k------>broj particija
%prvo probamo s tockama, jer slika izbacuje neku matricu s previse
%dimenzija :D

sigma=1;
X=rand(n/2,3)+10;
Y=rand(n/2,3)-10;

mat=[X;Y];

r=randperm(n);

for i = 1 : n
    PX(i,:) = mat(r(i),:);
end

PX
W = zeros( n, n );

for i = 1 : n
    for j = 1 : n
        norma  = norm( PX(i,:) - PX(j,:) );
        W( i, j ) = exp( - norma*norma/sigma );
    end
end

figure(2)
imagesc(PX);
%tu smo sad dobili tezinsku matricu W, pogledamo kako su rasporedjene

%do tu je ok, sad idemo na onaj algoritam

%1. korak
D=diag(W*ones(n,1));


%2. korak
%trazimo prvih K eigenvektora od (D^-1/2)W(D^-1/2), nakon toga racunamo Z
[V,S]=eig(sqrt(D)\W/(sqrt (D)));

[d,order] = sort(diag(S),'descend');  %# sort eigenvalues in descending order
V = V(:,order);

Z=sqrt(D)\V(:,1:k);

%3. korak
XZT=diag((diag(Z*Z')).^(-1/2))*Z;

%4. korak
i=randi(n,1);
R(:,1)=XZT(i,:)';    %5x1
c=zeros(n,1);     %100x1
for j=2:k
    c=c+abs(XZT*R(:,j-1));  %100x1
    [Y,I] = min(c);         
    R(:,j)=XZT(I,:)'
end

%5. korak
FIz=0;

while(1)
%6.korak
XT=XZT*R;

for i=1:n
    for l=1:k
        if (XT(i,l)==max(XT(i,1:k))) 
        XZ(i,l)=1;
        else 
        XZ(i,l)=0;
        end
    end
end

%7.korak
[U,O,UT]=svd(XZ'*XZT);
FIt=trace(O);
if abs(FIt-FIz)<1e-2
    break
end
abs(FIt-FIz)
FIz=FIt;
R=U*U';

end

figure(3)
imagesc(XZ);





