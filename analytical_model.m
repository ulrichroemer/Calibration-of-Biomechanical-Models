function [val] = analytical_model(theta,x)

     % theta lives in [0,1]
     theta(:,1) = 0.8 + 2*theta(:,1);
     theta(:,2) = 0.1 + 0.8*theta(:,2);
%      theta(:,3) = 1.5 + 2*theta(:,3);
%      theta(:,4) = 0 + 1.5*theta(:,4);

     val = zeros(length(x),size(theta,1));
     x = reshape(x,[length(x),1]);
     for i = 1:size(theta,1)
         val(:,i) = 2*theta(i,1)./theta(i,2).*x.^(theta(i,2));
         %+ 2*theta(i,3)./theta(i,4).*x.^(theta(i,4));     
     end
%     val = zeros(length(x),size(theta,1));
%     x = reshape(x,[length(x),1]);
%     for i = 1:size(theta,1)
%         val(:,i) = theta(i,2).*exp(theta(i,1)*x)-2;     
%     end
    
end