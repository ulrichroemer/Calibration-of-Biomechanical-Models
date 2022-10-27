function [response]= uqlab_KLE_surrogate(Vs,Es,mean_response,uqlab_surrogate_eval)

    response =  mean_response + Vs*Es.^(1/2)*uqlab_surrogate_eval';
    
end

