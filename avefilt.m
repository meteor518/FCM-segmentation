function d=avefilt(x,n)   
a(1:n,1:n)=1;   %a��n��nģ��,Ԫ��ȫ��1
p=size(x);   %����ͼ����p��q��,��p>n,q>n
x1=double(x);
x2=x1;
%A(a:b,c:d)��ʾA����ĵ�a��b��,��c��d�е�����Ԫ��
for i=1:p(1)-n+1
    for j=1:p(2)-n+1
        c=x1(i:i+(n-1),j:j+(n-1)).*a; %ȡ��x1�д�(i,j)��ʼ��n��n��Ԫ����ģ�����
        s=sum(sum(c));                 %��c����(��ģ��)�и�Ԫ��֮��
        x2(i+(n-1)/2,j+(n-1)/2)=s/(n*n); %��ģ���Ԫ�صľ�ֵ����ģ������λ�õ�Ԫ��
    end
end
d=uint8(x2);
