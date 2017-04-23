clc;clear all;tic
%% Input polynomial order
Order=19;

%% finding the coefficient of polynomial
InvCoe=her_POL(Order);
mainpoly=zeros(1,Order+1);
for i=1:Order+1
    mainpoly(i)=InvCoe(Order+2-i);
end
%% initial guess for bairstow algorithm
PolyOrder=Order;
Steps=Order-floor(Order/2)-1;                          %quadratic factors
even_recognizer=Order/2-floor(Order/2);
poly=mainpoly;
guess=findroots(Order);
first_aid_poly=zeros(1,Order+1);
second_aid_poly=zeros(1,Order+1);
r_array=ones(1,Steps)*(-0.5);
s_array=zeros(1,Steps);
rVector=zeros(1,Steps);
sVector=zeros(1,Steps);
for i=1:Steps
    s_array(Steps+1-i)=-guess(i)*guess(Order+1-i)+ 0.1;
end

%% bairstow method
for counter=1:Steps
    r(1)=r_array(counter);
    s(1)=s_array(counter);

    first_aid_poly(1)=1;
    first_aid_poly(2)=1;
    iter=1;
    while (abs(first_aid_poly(1))>10^(-12))||abs((first_aid_poly(2))>10^(-12))
    
        first_aid_poly(Order+1)=mainpoly(Order+1);
        first_aid_poly(Order)=mainpoly(Order)+r(iter)*first_aid_poly(Order+1);
        
        for i=Order-1:-1:1
            first_aid_poly(i)=mainpoly(i)+r(iter)*first_aid_poly(i+1)+s(iter)*first_aid_poly(i+2);
        end

        second_aid_poly(Order)=first_aid_poly(Order+1);
        second_aid_poly(Order-1)=first_aid_poly(Order)+r(iter)*second_aid_poly(Order);
        
        for i=Order-2:-1:1
            second_aid_poly(i)=first_aid_poly(i+1)+r(iter)*second_aid_poly(i+1)+s(iter)*second_aid_poly(i+2);
        end
        
        dr=(second_aid_poly(3)*first_aid_poly(1)-first_aid_poly(2)*second_aid_poly(2))/((second_aid_poly(2))^2-second_aid_poly(1)*second_aid_poly(3));
        ds=(first_aid_poly(2)*second_aid_poly(1)-first_aid_poly(1)*second_aid_poly(2))/((second_aid_poly(2))^2-second_aid_poly(1)*second_aid_poly(3));
        
        iter=iter+1;
        
        r(iter)=r(iter-1)+dr;
        s(iter)=s(iter-1)+ds;

        if iter>10000
            break
        end
        
    end
    
    rVector(counter)=r(iter);
    sVector(counter)=s(iter);
    
    for i=3:Order+1
        mainpoly(i-2)=first_aid_poly(i);
    end

    Order=Order-2;
    
end
%% finding the rooots 
for i=1:Steps
    delta=(rVector(i))^2+4*sVector(i);
    root(i,1)=(rVector(i)+sqrt(delta))/2;
    root(i,2)=(rVector(i)-sqrt(delta))/2;    
end

if even_recognizer==0
    delta2=mainpoly(2)^2-4*mainpoly(1)*mainpoly(3);
    root(Steps+1,1)=(-mainpoly(2)+sqrt(delta2))/(2*mainpoly(3));
    root(Steps+1,2)=(-mainpoly(2)-sqrt(delta2))/(2*mainpoly(3));     
else
    root(Steps+1,1)=-mainpoly(1)/mainpoly(2);
end

if (PolyOrder/2-floor(PolyOrder/2))==0
    root=sort(root(:));
else
    root=root(:);
    root(end)=[];
    root=sort(root);
end
%% showing the results
format long;
root=error_main(PolyOrder,root);
disp('-----------------------------  ' );
disp(['ROOTS OF HERMITE POLYNOMIAL ARE, H' num2str(PolyOrder) '(X)'] );
disp('-----------------------------  ' );
disp(root);

final_error=errors(PolyOrder,root);
disp('-----------------------------  ' );
disp( ['errors associated with each root respectively']);
disp('-----------------------------  ' );
disp([final_error'])
toc


