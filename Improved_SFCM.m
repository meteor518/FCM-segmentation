I=imread('r4.bmp');
subplot(311);imshow(I);
I=avefilt(I,3);%avefilt为均值滤波
I=double(I);
R=I(:,:,1);   %R,G,B为图像I的红绿蓝各分量
G=I(:,:,2);
B=I(:,:,3);
[m,n]=size(R);%m,n为R的行和列
c=5;    %c为聚类中心的个数
%*****统计Lab直方图的各像素数初始化*****
%将L,a,b分别以Lf,af,bf分档
Lf=10;
af=10;
bf=10;
L_num=100/Lf+1;
a_num=200/af+1;
b_num=200/bf+1;

for i=1:L_num                 
    for j=1:a_num
        for k=1:b_num
            HGLab(i,j,k)=0;%统计Lab直方图各点的像素数
        end
    end
end
%******************************************************************
%***********************将RGB分量Lab变换****************************
%******************************************************************
for i=1:m
    for j=1:n
        W=[0.433910 0.376220 0.189860;0.212649 0.715169 0.072182;0.017756 0.109478 0.872915]*[R(i,j)/255;G(i,j)/255;B(i,j)/255];%变换的固定公式,得到3*1的矩阵
        X=W(1,1);
        Y=W(2,1);               %X,Y,Z对应取W矩阵的每一行的值,方便以下Lab变换使用
        Z=W(3,1);               
        if(Y>0.0088556)         %Lab变换后的 L值
            L(i,j)=(116*Y.^(1/3)-16); 
        else
            L(i,j)=903.3*Y;
        end
        a(i,j)=500.0*(f(X)-f(Y));      %Lab变换后的 a值
        b(i,j)=200.0*(f(Y)-f(Z));      %Lab变换后的 b值
    end
end

%******************************************************************
%***********************直方图统计并找峰值**************************
%******************************************************************
%直方图统计
for i=1:m
    for j=1:n
        HGLab(fix(L(i,j)/Lf)+1,fix((a(i,j)+100)/af)+1,fix((b(i,j)+100)/bf)+1)=HGLab(fix(L(i,j)/Lf)+1,fix((a(i,j)+100)/af)+1,fix((b(i,j)+100)/bf)+1)+1;
    end
end
%找峰值
for i=1:L_num
    for p=1:2
        Labmax(i,p)=HGLab(i,1,1);
        L_max(i,p)=i;  
        a_max(i,p)=1;   
        b_max(i,p)=1;   
        for j=1:a_num
           for k=1:b_num
               if(HGLab(i,j,k)>=Labmax(i,p))
                   Labmax(i,p)=HGLab(i,j,k);
                   a_max(i,p)=j;
                   b_max(i,p)=k;
               end
           end
        end
        HGLab(i,a_max(i,p),b_max(i,p))=0;
    end
end        
%找出峰值中前c个最大值
ii=1;
while(ii<c+1)
    Lab(ii,1)=Labmax(1,1);
    for i=1:L_num  
        for p=1:2
            if(Lab(ii,1)<=Labmax(i,p))
                Lab(ii,1)=Labmax(i,p);
                Lmax(ii,1)=L_max(i,p);
                amax(ii,1)=a_max(i,p);
                bmax(ii,1)=b_max(i,p);
                i_max=i;
                p_max=p;
            end
        end
    end
    Labmax(i_max,p_max)=0;
    ii=ii+1;
end
%将最大值对应的Lab分档变换为初始聚类中心
z=1;
for i=1:ii-1
    LL=(Lmax(i,1)-0.5)*Lf;
    aa=(amax(i,1)-0.5-1/af)*af;            
    bb=(bmax(i,1)-0.5-1/bf)*bf;

    Choose_L(z,1)=LL;
    Choose_a(z,1)=aa;
    Choose_b(z,1)=bb;

    z=z+1;
end

%******************************************************************
%**********************在Lab空间下进行SFCM算法**********************
%******************************************************************
e=1;        %停止阈值
count=0;    %迭代次数
flag=true;  %迭代循环标志
A=0.32;      %系数因子
sigma=50;    %标准差

SumWL=0;
SumWa=0;
SumWb=0;
SumWLL=0;
SumWaa=0;
SumWbb=0;
%**********求权重WLab**********
for i=2:m-1
    for j=2:n-1            
        for t1=1:3
            for t2=1:3
                WL(t1,t2)=exp(-(L(i,j)-L(i-2+t1,j-2+t2))^2/(2.0*sigma^2));
                Wa(t1,t2)=exp(-(a(i,j)-a(i-2+t1,j-2+t2))^2/(2.0*sigma^2));
                Wb(t1,t2)=exp(-(b(i,j)-b(i-2+t1,j-2+t2))^2/(2.0*sigma^2));
            end
        end
        for t1=1:3
            for t2=1:3
                SumWL=SumWL+WL(t1,t2);
                SumWa=SumWa+Wa(t1,t2);
                SumWb=SumWb+Wb(t1,t2);

                SumWLL=SumWLL+WL(t1,t2)*L(i,j);
                SumWaa=SumWaa+Wa(t1,t2)*a(i,j);
                SumWbb=SumWbb+Wb(t1,t2)*b(i,j);
            end
        end
        L1(i,j)=SumWLL/SumWL;
        a1(i,j)=SumWaa/SumWa;
        b1(i,j)=SumWbb/SumWb;
    end
end
while(flag&&count<1000)   
    for i=2:m-1
        for j=2:n-1 
            tL(i,j)=(L(i,j)+A*L1(i,j))/(1+A);
            ta(i,j)=(a(i,j)+A*a1(i,j))/(1+A);
            tb(i,j)=(b(i,j)+A*b1(i,j))/(1+A);
            
            %**********求距离dist**********
            for k=1:c
                distL(k)=(Choose_L(k,1)-tL(i,j))*(Choose_L(k,1)-tL(i,j));
                dista(k)=(Choose_a(k,1)-ta(i,j))*(Choose_a(k,1)-ta(i,j));
                distb(k)=(Choose_b(k,1)-tb(i,j))*(Choose_b(k,1)-tb(i,j));

                distLab(k)=(distL(k)+dista(k)+distb(k))/3.0;
            end

            %**********求隶属度LabU**********
            p=(i-1)*n+j;%每个像素点的位置，p表示是第几个像素点

            %先对距离进行判别，若有0距离，则该U为1，其他为0；若没有按公式求
            for k=1:c
                if(distLab(k)==0)
                    LabU(k,p)=1;
                    if(k==c)
                        k=0;
                    end
                    break;
                end
            end
            if(k<c)     %代表上述有距离为0的点
                for q=1:c
                    if(k==0)
                        if(q>c||q<c)
                            LabU(q,p)=0;
                        end
                    else
                        if(q>k||q<k)
                            LabU(q,p)=0;
                        end
                    end
                end
            else        %没有距离为0的点
                Sum_dist=0;
                for r=1:c
                    Sum_dist=Sum_dist+(1/distLab(r));
                end
                for k=1:c
                    LabU(k,p)=(1/distLab(k))*(1/Sum_dist);
                end
            end
        end
    end

    %*********求变换后的聚类中心矩阵V**********
    for k=1:c
        VL=0;
        Va=0;
        Vb=0;
        V=0;
        for i=2:m-1
            for j=2:n-1
                p=(i-1)*n+j;%每个像素点的位置，p表示是第几个像素点
                VL=VL+LabU(k,p)*LabU(k,p)*tL(i,j);
                Va=Va+LabU(k,p)*LabU(k,p)*ta(i,j);
                Vb=Vb+LabU(k,p)*LabU(k,p)*tb(i,j);
                V =V +LabU(k,p)*LabU(k,p);
            end
        end
        Old_L(k,1)=Choose_L(k,1);
        Old_a(k,1)=Choose_a(k,1);   %先将原来的聚类中心点保存在Old中
        Old_b(k,1)=Choose_b(k,1);

        Choose_L(k,1)=VL/V;
        Choose_a(k,1)=Va/V;         %新的聚类中心点
        Choose_b(k,1)=Vb/V;
    end
    flag=false;
    for k=1:c
        if(abs(Choose_L(k,1)-Old_L(k,1))>e||abs(Choose_a(k,1)-Old_a(k,1))>e||abs(Choose_b(k,1)-Old_b(k,1))>e)
            count=count+1;
            flag=true;
            break;
        end
    end
end

%***************判断C个中心中的焊盘色****************
for k=1:c
    if(Choose_L(k,1)>55&&Choose_L(k,1)<85 &&Choose_a(k,1)>-11&&Choose_a(k,1)<11&&Choose_b(k,1)>-11&&Choose_b(k,1)<11)
        gray_num=k;
    end
end
    
%*********************将Choose_L,a,b变换为R,G,B*********************
%*********************将Choose_L,a,b变换为R,G,B*********************
for k=1:c
    fy = (Choose_L(k,1)+16.0)/116.0;
    fy = fy*fy*fy;

    if(fy > 0.008856)
        y=fy;
    else
        fy = Choose_L(k,1)/903.3;
    end

    if(fy > 0.008856)
        fy = fy.^(1.0/3.0);
    else
        fy = 7.787*fy+16.0/116.0;
    end

    fx = Choose_a(k,1)/500.0 + fy;
    if(fx > 0.206893)
        x = fx.^3.0;
    else
        x = (fx-16.0/116.0)/7.787;
    end

    fz = fy - Choose_b(k,1)/200.0;
    if(fz > 0.206893)
        z = fz.^3;
    else
        z = (fz-16.0/116.0)/7.787;
    end

    x = x*0.950456*255.0;
    y = y*255.0;
    z = z*1.088754*255.0;

    r_rgb(k,1) =  3.240479*x  - 1.537150*y - 0.498535*z;
    g_rgb(k,1)=  -0.969256*x + 1.875992*y + 0.041556*z;
    b_rgb(k,1)=  0.055648*x  - 0.204043*y + 1.057311*z;
end
%*********分割图像显示，留下焊盘色，其他为白色*************
for i=1:m
    for j=1:n
        for k=1:c
            rr(k)=abs(R(i,j)-r_rgb(k,:));
            gg(k)=abs(G(i,j)-g_rgb(k,:));
            bb(k)=abs(B(i,j)-b_rgb(k,:));
            e1(k)=(rr(k)+gg(k)+bb(k))/3.0;
            rg(k)=abs(abs(r_rgb(k,:)-g_rgb(k,:))-abs(R(i,j)-G(i,j)));
            rb(k)=abs(abs(r_rgb(k,:)-b_rgb(k,:))-abs(R(i,j)-B(i,j)));
            gb(k)=abs(abs(g_rgb(k,:)-b_rgb(k,:))-abs(G(i,j)-B(i,j)));
            e2(k)=(rg(k)+rb(k)+gb(k))/3.0;
            if(k==gray_num)
                e(k)=(e1(k)+e2(k))/2.0;
            else
                wrg(k)=abs(r_rgb(k,:)-g_rgb(k,:))/(abs(r_rgb(k,:)-g_rgb(k,:))+abs(r_rgb(k,:)-b_rgb(k,:))+abs(b_rgb(k,:)-g_rgb(k,:)));
                wrb(k)=abs(r_rgb(k,:)-b_rgb(k,:))/(abs(r_rgb(k,:)-g_rgb(k,:))+abs(r_rgb(k,:)-b_rgb(k,:))+abs(b_rgb(k,:)-g_rgb(k,:)));
                wgb(k)=abs(b_rgb(k,:)-g_rgb(k,:))/(abs(r_rgb(k,:)-g_rgb(k,:))+abs(r_rgb(k,:)-b_rgb(k,:))+abs(b_rgb(k,:)-g_rgb(k,:)));
                e3(k)=wrg(k)*rg(k)+wrb(k)*rb(k)+wgb(k)*gb(k);
                e(k)=(e1(k)+e3(k))/2.0;    
            end
        end
        p=1;
        for k=1:c-1
            if(e(p)>e(k+1))
                p=k+1;
            else
                p=p;
            end
        end
        R(i,j)=r_rgb(p,:);
        G(i,j)=g_rgb(p,:);
        B(i,j)=b_rgb(p,:);
        if(p>gray_num||p<gray_num)
            R(i,j)=255;
            G(i,j)=255;
            B(i,j)=255;
        else
            R(i,j)=0;
            G(i,j)=0;
            B(i,j)=0;
        end
        I(i,j,1)=R(i,j);
        I(i,j,2)=G(i,j);
        I(i,j,3)=B(i,j);
    end
end
I=uint8(I);
subplot(312);imshow(I);

%********计算目标函数Jm的值***********
Jm=0;
for k=1:c
    for i=2:m-1
        for j=2:n-1
            p=(i-1)*n+j;
            dL=(Choose_L(k,1)-tL(i,j))*(Choose_L(k,1)-tL(i,j));
            da=(Choose_a(k,1)-ta(i,j))*(Choose_a(k,1)-ta(i,j));
            db=(Choose_b(k,1)-tb(i,j))*(Choose_b(k,1)-tb(i,j));

            dLab=(dL+da+db)/3.0;
            Jm=Jm+LabU(k,p)*LabU(k,p)*dLab;
        end
    end
end 

%********比较新的聚类中心与原始的距离***********
for k=1:c
    L1=(Choose_L(k,1)-Old_L(k,1))*(Choose_L(k,1)-Old_L(k,1));
    a1=(Choose_a(k,1)-Old_a(k,1))*(Choose_a(k,1)-Old_a(k,1));
    b1=(Choose_b(k,1)-Old_b(k,1))*(Choose_b(k,1)-Old_b(k,1));
    d(k,1)=(L1+a1+b1)/3;
end
D=0;
for k=1:c
    D=D+d(k,1);
end
D=D/c;

Area1=0;
Area2=0;     %焊盘面积
for i=1:m
    if(i>=85&&i<=157)
        Width(i,1)=0;
        Width(i,2)=0;
        for j=1:n
            if(j>=80&&j<=130)
                if(R(i,j)==0)
                    R(i,j)=255;
                    I(i,j,1)=R(i,j);
                    I(i,j,2)=G(i,j);
                    I(i,j,3)=B(i,j);
                    Area1=Area1+1;
                    Width(i,1)=Width(i,1)+1;
                end
            end
            if(j>=260&&j<=310)
                if(R(i,j)==0)
                    R(i,j)=255;
                    I(i,j,1)=R(i,j);
                    I(i,j,2)=G(i,j);
                    I(i,j,3)=B(i,j);
                    Area2=Area2+1;
                    Width(i,2)=Width(i,2)+1;
                end
            end
        end             
    end
end
I=uint8(I);
subplot(313);imshow(I);

Width1=Width(85,1);
Width2=Width(85,2);  %焊盘宽度
for i=85:157
    if(Width1<Width(i,1))
        Width1=Width(i,1);
    end
    if(Width2<Width(i,2))
        Width2=Width(i,2);
    end
end

for j=1:n
    if(j>=80&&j<=130||j>=260&&j<=310)
        Height(1,j)=0;
        for i=1:m
            if(i>=85&&i<=157)
                if(G(i,j)==0)
                    Height(1,j)=Height(1,j)+1;
                end
            end
        end
    end
end
Height1=Height(1,80);
Height2=Height(1,260);  %焊盘高度
for j=80:130
    if(Height1<Height(1,j))
        Height1=Height(1,j);
    end
end
for j=260:310
    if(Height2<Height(1,j))
        Height2=Height(1,j);
    end
end