 N=10
for i = 1:N       
     p(i) = i;
 end
p
for i = 1:N       
     %r = i + (rand() % (52-i));
     
     r=i+randi([0,N-i],1,1);
     
     temp = p(i); 
     p(i) = p(r); 
     p(r) = temp;
 end
p