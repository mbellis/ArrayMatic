n1=[126:140];
n2=[196:210];
n3=[212:226];
h1=figure;
h2=figure;
h3=figure;
for i=1:15
    net=sprintf('n%u',n1(1)+i-1);
    eval(sprintf('cd /home/mbellis/array1/sosma/net/m8/%s',net))
    a1=load_data(sprintf('a_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    c1=load_data(sprintf('c_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    f1=load_data(sprintf('f_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);


    net=sprintf('n%u',n2(1)+i-1);
    eval(sprintf('cd /home/mbellis/array1/sosma/net/m8/%s',net))
    a2=load_data(sprintf('a_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    c2=load_data(sprintf('c_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    f2=load_data(sprintf('f_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);

    net=sprintf('n%u',n3(1)+i-1);
    eval(sprintf('cd /home/mbellis/array1/sosma/net/m8/%s',net))
    a3=load_data(sprintf('a_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    c3=load_data(sprintf('c_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    f3=load_data(sprintf('f_m8_%s.4mat',net),'./',45101,45101,'uint8','ieee-le',1:45101,1000);
    figure(h1)
    subplot(3,5,i)
    hold on
    plot(a1,a3,'b.')
    plot(a1,a2,'c.')
    title(sprintf('n%u n%u n%u',n1(1)+i-1,n2(1)+i-1,n3(1)+i-1))

    figure(h2)
    subplot(3,5,i)
    hold on
    plot(c1,c3,'r.')
    plot(c1,c2,'m.')
    title(sprintf('n%u n%u n%u',n1(1)+i-1,n2(1)+i-1,n3(1)+i-1))

    figure(h3)
    subplot(3,5,i)
    hold on
    plot(f1,f2,'k.')
    plot(f1,f3,'g.')
    title(sprintf('n%u n%u n%u',n1(1)+i-1,n2(1)+i-1,n3(1)+i-1))

end