function [ConvA Convb ConcA Concb] = TrilinearEnvelopeEngine(bnd, wIndex,Dim)
%output: The trilinear convex and concave envelope for w=x1x2x3
%        ConvA*x >= Convb
%        ConcA*x >= Concb
%Input: the bounds on x1,x2,x3. bound = [x1_L x1_U Index1
%                                        x2_L x2_U Index2
%                                        x3_L x3_U Index3];
%       wIndex: the variable index of w
%       Dim: Total number of variabls
% Version: 2011.08.03

%The number of Intervals across 0
NumCross = 0;
NumPos = 0;
NumNeg = 0;
for i = 1:3
    if bnd(i,1)*bnd(i,2) < 0
        NumCross = NumCross +1;
    end
end

for i = 1:3
    if bnd(i,2) <= 0
        NumNeg = NumNeg + 1;
    end
end

for i = 1:3
    if bnd(i,1) >= 0
        NumPos = NumPos + 1;
    end
end

if NumCross + NumPos + NumNeg ~= 3
    fprintf('\n\n******** There must be an error in trilinear envelope********\n\n');
    pause;
end

label_convex = 'false';
counter_convex = 0;
label_concave = 'true';
counter_concave = 0;

switch NumCross
    %Positive or Negative domains
    case 0
        switch NumPos
            % three positive domains
            case 3
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
                       (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) && ...
                       (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
                       (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                   
                        %convex envelope: there are 6 constraints 
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(2) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);    
                        
                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(3) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            
                        ConvA(5,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                        ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(5) = -theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                               
                        theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(6,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                        ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(6) = -theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        %concave envelope: there are 6 constraints
                        ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(1) =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(2) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                  +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(3) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                              
                        ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(4) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                              
                        ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(5) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                        ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(6) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                              
                        label_convex = 'true';
                        counter_convex = counter_convex + 1;
                        label_concave = 'true';
                        counter_concave = counter_concave + 1;
                        
                        break;                        
                    end
                end
                            
            %two positive domains, one negative domain
            case 2
               
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if  bnd(pmt(3,i),2) <= 0  && bnd(pmt(1,i),1) >= 0 && bnd(pmt(2,i),1) >= 0
                        
                        ConvA = sparse(zeros(6,Dim));
                        Convb = sparse(zeros(6,1));
                        %convex envelope: there are 6 constraints 
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(1) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...;
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                               
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(2) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);    
                        
                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(3) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        
                                                  
                        ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(5) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                               
                        
                        ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(6) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                               
                        %concave envelope         
                        if (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                           (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) && ...
                           (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                           (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                            %%%%
                            ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                            ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    
                        
                            ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Concb(3) =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                            ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Concb(4) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            
                            ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConcA(5,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2));
                            Concb(5) =  theta*bnd(pmt(3,i),2)/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2))...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                            Concb(6) =  theta*bnd(pmt(3,i),1)/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1))...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            
                                                          
                            label_convex = 'true';
                            counter_convex = counter_convex + 1;
                            label_concave = 'true';
                            counter_concave = counter_concave + 1;
                        
                            break;
                               
                        elseif (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                               (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                            %%%%
                            ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                            ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    
                        
                            ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Concb(3) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                            ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Concb(4) =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        
                            theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            
                            ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(5,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                            ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Concb(5) =  theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                            ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            Concb(6) =  theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);  
                                                          
                            label_convex = 'true';
                            counter_convex = counter_convex + 1;
                            label_concave = 'true';
                            counter_concave = counter_concave + 1;
                        
                            break;
                        else
                            continue;
                        end
                    end
                end                
            
            % one positive domain, two negative domain
            case 1
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
            
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if  bnd(pmt(1,i),1) >= 0 && bnd(pmt(2,i),2) <= 0 && bnd(pmt(3,i),2) <= 0
                        
                        ConcA = sparse(zeros(6,Dim));
                        Concb = sparse(zeros(6,1));
                        %concave envelope
                        ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(1) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(2) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                  +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(3) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                              
                        ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(4) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                              
                        ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(5) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(6) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                              
                        %convex envelope: there are 6 constraints 
                        if (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
                           (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) && ...
                           (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
                           (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                        
                            ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Convb(2) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);    

                            ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(3) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            ConvA(5,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                            ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(5) = -theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConvA(6,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                            ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(6) = -theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                                           
                            label_convex = 'true';
                            counter_convex = counter_convex + 1;
                            label_concave = 'true';
                            counter_concave = counter_concave + 1;
                        
                             break;
                             
                        elseif (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) <= ...
                               (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) && ...
                               (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) <= ...
                               (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                            ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Convb(2) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);    

                            ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(3) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(4) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConvA(5,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                            ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(5) = -theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                            ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(6) = -theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                                                            
                            label_convex = 'true';
                            counter_convex = counter_convex + 1;
                            label_concave = 'true';
                            counter_concave = counter_concave + 1;
                            
                            break;
                        else
                           continue;
                        end
                    end                               
                end                
            
            %three negative domains
            case 0
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if   (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                         (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) && ...
                         (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                         (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))      
                     
                        %convex envelope: there are 6 constraints 
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(1) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...;
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                               
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(2) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    
                        
                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(3) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(4) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                                                  
                        ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(5) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                               
                        
                        ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(6) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                               
                        
                        %concave envelope         
                        ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);    
                        
                        ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(3) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                        ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(4) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            
                        ConcA(5,bnd(pmt(1,i),3)) = theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                        ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(5) =  theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                               
                        theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(1,i),3)) = theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                        ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(6) =  theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                                                         
                        label_convex = 'true';
                        counter_convex = counter_convex + 1;
                        label_concave = 'true';
                        counter_concave = counter_concave + 1;
                        
                        break;                           
                    end
                end
        end                 

    %mixed domains: 1 cross    
    case 1
        switch NumPos
            % one cross-zero domain, two positive domains
            case 2 %verified on 2011.08.03
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(3,i),1) <= 0 &&  bnd(pmt(3,i),2) >= 0 && bnd(pmt(1,i),1) >= 0 && bnd(pmt(2,i),1) >= 0 %corrected 2011.08.03
                        %convex envelope
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(1) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                                   
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(2) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);    

                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(3) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                        ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(5) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                        theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(3,i),3)) = -theta/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                        Convb(6) = -theta*bnd(pmt(3,i),1)/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1))...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                               
                        %concave envelope
                        ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(1) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(2) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);    
                        
                        ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(3) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                        ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(4) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
   
                        ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(5) =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                               
                        theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(6,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2));
                        Concb(6) =  theta*bnd(pmt(3,i),2)/(bnd(pmt(3,i),1)-bnd(pmt(1,i),2))...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                                                             
                        label_convex = 'true';
                        counter_convex = counter_convex + 1;
                        label_concave = 'true';
                        counter_concave = counter_concave + 1;
                        
                        break;
                    end
                end

            % one cross-zero domain, one positive domain and one negative domain 
            case 1 
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(1,i),1) >= 0 &&  bnd(pmt(3,i),2) <= 0 && bnd(pmt(2,i),1) <= 0 && bnd(pmt(2,i),2) >= 0
                        %convex envelope
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(1) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                   
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(2) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);    

                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(3) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);


                        ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(5) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);


                        theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                        ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(6) = -theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                               
                        %concave envelope
                        ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(1) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        
                        ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(2) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);    
                        
                        ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(3) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(4) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
   
                        ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(5) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                               
                        theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                        ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(6) =  theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                                             
                        label_convex = 'true';
                        counter_convex = counter_convex + 1;
                        label_concave = 'true';
                        counter_concave = counter_concave + 1;
                        
                        break;
                    end
                end
                
            % one cross-zero domain, two negative domains
            case 0
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(2,i),2) <= 0 &&  bnd(pmt(3,i),2) <= 0 && bnd(pmt(1,i),1) <= 0 && bnd(pmt(1,i),2) >= 0
                        %convex envelope
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(1) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                   
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(2) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    

                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(3) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);


                        ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(5) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);


                        theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                        ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(6) = -theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                               
                        %concave envelope
                        ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        
                        ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(2) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);    
                        
                        ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Concb(3) =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        
                        ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Concb(4) =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
   
                        ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(5) =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                               
                        theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(1,i),3)) = theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                        ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Concb(6) =  theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                                             
                        label_convex = 'true';
                        counter_convex = counter_convex + 1;
                        label_concave = 'true';
                        counter_concave = counter_concave + 1;
                        
                        break;
                    end
                end                
        end
        
    %mixed domains: 2 cross    
    case 2
        switch NumPos
            % two cross-zero domains, one positive domain
            case 1  %checked on 2011.08.03
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(1,i),1) >= 0 && bnd(pmt(2,i),1) <= 0 && bnd(pmt(2,i),2) >= 0 && bnd(pmt(3,i),1) <= 0 && bnd(pmt(3,i),2) >= 0 &&...
                       (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) <= ...
                       (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                   
                        %convex envelope
                        ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                        ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                        Convb(1) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                                   
                        ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                        ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                        Convb(2) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);    

                        ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                        Convb(3) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                        ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                        ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                        ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(4) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                        theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(5,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                        ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                        Convb(5) = -theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);


                        theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                        ConvA(6,bnd(pmt(3,i),3)) = -theta/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                        Convb(6) = -theta*bnd(pmt(3,i),1)/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1))...
                                   -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                   -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                   +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                                                             
                        label_convex = 'true';
                        counter_convex = counter_convex + 1;

                        break;
                    end
                end
                
                %concave envelope
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(1,i),1) >= 0 && bnd(pmt(2,i),1) <= 0 && bnd(pmt(2,i),2) >= 0 && bnd(pmt(3,i),1) <= 0 && bnd(pmt(3,i),2) >= 0
                       if (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                          (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))
                   
                        
                            ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Concb(1) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    

                            ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Concb(3) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                      +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Concb(4) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                      +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConcA(5,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2));
                            Concb(5) = theta*bnd(pmt(3,i),2)/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2))...
                                      +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                      +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                      -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);


                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(6,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                            ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Concb(6) =  theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                                   
                            label_concave = 'true';
                            counter_concave = counter_concave + 1;
                            
                            break;
                       else
                            ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Concb(1) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    

                            ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Concb(3) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                      +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Concb(4) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                      +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConcA(5,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                            Concb(5) = theta*bnd(pmt(3,i),1)/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1))...
                                      +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                      +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                      -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConcA(6,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                            ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Concb(6) =  theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                   
                            label_concave = 'true';
                            counter_concave = counter_concave + 1;
                            
                            break;
                       end
                    end
                end   
   
            % two cross-zero domains, one negative domain
            case 0 
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(3,i),2) <= 0 && bnd(pmt(1,i),1) <= 0 && bnd(pmt(1,i),2) >= 0 && bnd(pmt(2,i),1) <= 0 && bnd(pmt(2,i),2) >= 0
                        if (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
                           (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                   
                            %convex envelope
                            ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Convb(2) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);    

                            ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(3) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(4) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(5,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                            ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(5) = -theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                            theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                            ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(6) = -theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                                   
                            label_convex = 'true';
                            counter_convex = counter_convex + 1;
                        
                            break;
                        else
                            ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                            ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                            Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                            ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                            ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                            Convb(2) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);    

                            ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                            ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(3) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                            Convb(4) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                            theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                            ConvA(5,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                            ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                            ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(5) = -theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);


                            theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                    -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                    +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                            ConvA(6,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                            ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                            Convb(6) = -theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                                       -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                                       -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                       +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                            label_convex = 'true';
                            counter_convex = counter_convex + 1;
                            break;
                        end
                    end
                end
                
                %concave envelope
                pmt = [1 1 2 2 3 3
                       2 3 1 3 2 1
                       3 2 3 1 1 2];
                for i = 1:6
                    if bnd(pmt(3,i),2) <= 0 && bnd(pmt(1,i),1) <= 0 && bnd(pmt(1,i),2) >= 0 && bnd(pmt(2,i),1) <= 0 && bnd(pmt(2,i),2) >= 0 && ...
                       (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                       (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                   
                        
                       ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                       ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                       ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                       Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                       ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                       ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                       ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                       Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    

                       ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                       ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                       ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                       Concb(3) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                 +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                       ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                       ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                       ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                       Concb(4) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                 +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                       theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                               -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                               -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                               +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                       ConcA(5,bnd(pmt(1,i),3)) = theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                       ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                       ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                       Concb(5) = theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                                 +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                 +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                                 -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                       theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                               -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                               -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                               +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                       ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                       ConcA(6,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                       ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                       Concb(6) =  theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                                  +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                                  +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                                  -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                       label_concave = 'true';
                      counter_concave = counter_concave + 1;
                      break;
                    end
                end
        end
    %mixed domains: 3 cross    
    case 3 %verified on 2011.08.03 
        %%%%%
        %convex envelope        
        pmt = [1 1 2 2 3 3
               2 3 1 3 2 1
               3 2 3 1 1 2];
        for i = 1:6
            if (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))&&...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) <= ...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2) + 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1))&&...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
               (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2) + 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1))&&...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
               (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))
                %convex envelope
                ConvA = sparse(zeros(5,Dim));
                Convb = sparse(zeros(5,1));
                
                ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Convb(2) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);    

                ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Convb(3) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Convb(4) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
          
                thetax = 1/2*(bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)-...
                              bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)-bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))/...
                              (bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                thetay = 1/2*(bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)-...
                              bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)-bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))/...
                              (bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                thetaz = 1/2*(bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)-...
                              bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)-bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))/...
                              (bnd(pmt(3,i),1)-bnd(pmt(3,i),2));
                thetac = bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) - thetax*bnd(pmt(1,i),1)-thetay*bnd(pmt(2,i),1)-thetaz*bnd(pmt(3,i),1);
                
                ConvA(5,bnd(pmt(1,i),3)) = -thetax;
                ConvA(5,bnd(pmt(2,i),3)) = -thetay;
                ConvA(5,bnd(pmt(3,i),3)) = -thetaz;
                Convb(5) = thetac;
                
                label_convex = 'true';
                counter_convex = counter_convex + 1;
                            
                break;
                
            elseif (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
                   (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))
                %convex envelope
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Convb(2) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);    

                ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Convb(3) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                
                theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);                         
                ConvA(4,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                ConvA(4,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Convb(4) = -theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConvA(5,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConvA(5,bnd(pmt(3,i),3)) = -theta/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                Convb(5) = -theta*bnd(pmt(3,i),1)/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1))...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConvA(6,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                ConvA(6,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Convb(6) = -theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
               label_convex = 'true';
               counter_convex = counter_convex + 1;
                break;
                
            elseif (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) >= ...
                   (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2) + 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1))
                %convex envelope
                ConvA = sparse(zeros(6,Dim));
                Convb = sparse(zeros(6,1));
                
                ConvA(1,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConvA(1,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConvA(1,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Convb(1) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                ConvA(2,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConvA(2,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConvA(2,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Convb(2) = -2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);    

                ConvA(3,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConvA(3,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConvA(3,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Convb(3) = -2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                
                theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);                         
                ConvA(4,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConvA(4,bnd(pmt(2,i),3)) = -theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                ConvA(4,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Convb(4) = -theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))... %corrected here on 2011.08.03
                           -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConvA(5,bnd(pmt(1,i),3)) = -theta/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                ConvA(5,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConvA(5,bnd(pmt(3,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Convb(5) = -theta*bnd(pmt(1,i),1)/(bnd(pmt(1,i),2)-bnd(pmt(1,i),1))...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);


                theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConvA(6,bnd(pmt(1,i),3)) = -bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConvA(6,bnd(pmt(2,i),3)) = -bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConvA(6,bnd(pmt(3,i),3)) = -theta/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2));
                Convb(6) = -theta*bnd(pmt(3,i),2)/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2))...
                           -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);   
                label_convex = 'true';
                counter_convex = counter_convex + 1;
               break;
            else
                continue;
            end
        end
        
        %concave envelope        
        pmt = [1 1 2 2 3 3
               2 3 1 3 2 1
               3 2 3 1 1 2];
        for i = 1:6
            if (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) >= ...
               (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2) + 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1))&&...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
               (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))&&...
               (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1) + 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2))&&...
               (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) >= ...
               (bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2) + 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1))
                %concave envelope
                ConcA = sparse(zeros(5,Dim));
                Concb = sparse(zeros(5,1));
                
                ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    

                ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Concb(3) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);

                ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Concb(4) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
          
                thetax = 1/2*(bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)-...
                              bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)-bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1))/...
                              (bnd(pmt(1,i),2)-bnd(pmt(1,i),1));
                thetay = 1/2*(bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)-...
                              bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)-bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1))/...
                              (bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                thetaz = 1/2*(bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)-...
                              bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)-bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1))/...
                              (bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                thetac = bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2) - thetax*bnd(pmt(1,i),2)-thetay*bnd(pmt(2,i),2)-thetaz*bnd(pmt(3,i),2);
                
                ConcA(5,bnd(pmt(1,i),3)) = thetax;
                ConcA(5,bnd(pmt(2,i),3)) = thetay;
                ConcA(5,bnd(pmt(3,i),3)) = thetaz;
                Concb(5) = -thetac;
                
                label_concave = 'true';
                counter_concave = counter_concave + 1;
                
                break;
                
            elseif (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)) <= ...
                   (bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2) + 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1))
                %convex envelope
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Concb(1) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);

                ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);    

                ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Concb(3) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                
                theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);                         
                ConcA(4,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConcA(4,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2));
                Concb(4) =  theta*bnd(pmt(3,i),2)/(bnd(pmt(3,i),1)-bnd(pmt(3,i),2))...
                           +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);

                theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConcA(5,bnd(pmt(1,i),3)) = theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                ConcA(5,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Concb(5) =  theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);


                theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConcA(6,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2));
                ConcA(6,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Concb(6) =  theta*bnd(pmt(2,i),2)/(bnd(pmt(2,i),1)-bnd(pmt(2,i),2))...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                
                label_concave = 'true';
                counter_concave = counter_concave + 1;
                
                break;
                
            elseif (bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)+bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)+bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)) <= ...
                   (bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1) + 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2))
                %concave envelope
                ConcA = sparse(zeros(6,Dim));
                Concb = sparse(zeros(6,1));
                
                ConcA(1,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                ConcA(1,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),1);
                ConcA(1,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),1);
                Concb(1) = 2*bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                ConcA(2,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConcA(2,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),1);
                ConcA(2,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),2);
                Concb(2) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);    

                ConcA(3,bnd(pmt(1,i),3)) = bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                ConcA(3,bnd(pmt(2,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(3,i),2);
                ConcA(3,bnd(pmt(3,i),3)) = bnd(pmt(1,i),2)*bnd(pmt(2,i),1);
                Concb(3) = 2*bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2);
                
                theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);                         
                ConcA(4,bnd(pmt(1,i),3)) = theta/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2));
                ConcA(4,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConcA(4,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Concb(4) =  theta*bnd(pmt(1,i),2)/(bnd(pmt(1,i),1)-bnd(pmt(1,i),2))...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);

                theta =  bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1);
                ConcA(5,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConcA(5,bnd(pmt(2,i),3)) = theta/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1));
                ConcA(5,bnd(pmt(3,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(2,i),2);
                Concb(5) =  theta*bnd(pmt(2,i),1)/(bnd(pmt(2,i),2)-bnd(pmt(2,i),1))...
                           +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),1)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);


                theta =  bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                        -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1)...
                        -bnd(pmt(1,i),1)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                        +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConcA(6,bnd(pmt(1,i),3)) = bnd(pmt(2,i),2)*bnd(pmt(3,i),2);
                ConcA(6,bnd(pmt(2,i),3)) = bnd(pmt(1,i),1)*bnd(pmt(3,i),2);
                ConcA(6,bnd(pmt(3,i),3)) = theta/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1));
                Concb(6) =  theta*bnd(pmt(3,i),1)/(bnd(pmt(3,i),2)-bnd(pmt(3,i),1))...
                           +bnd(pmt(1,i),1)*bnd(pmt(2,i),1)*bnd(pmt(3,i),2)...
                           +bnd(pmt(1,i),2)*bnd(pmt(2,i),2)*bnd(pmt(3,i),2)...
                           -bnd(pmt(1,i),2)*bnd(pmt(2,i),1)*bnd(pmt(3,i),1);
                label_concave = 'true';
                counter_concave = counter_concave + 1;
               break;
            else
                continue;
            end
        end      
end

if strcmp(label_convex,'false') || strcmp(label_concave,'false')|| counter_convex ~= 1 || counter_concave ~= 1 
    fprintf('\n\n******** There must be an error in trilinear envelope********\n\n');
    pause;
end

% set the coefficient of variable w to be 1 or -1
ConvA(:,wIndex) = ones(size(ConvA,1),1);
ConcA(:,wIndex) = -ones(size(ConcA,1),1);
end
 
        

