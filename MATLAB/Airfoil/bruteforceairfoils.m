%Generates array of airfoil numbers to check GHstep




function Airfoils = bruteforceairfoils(A,B,CC)

    Airfoils = cell(2,4);
    
    Airfoils{1,1} = [num2str(A) num2str(B) num2str(CC)];
    Airfoils{1,2} = [num2str(A+1) num2str(B) num2str(CC)];
    Airfoils{1,3} = [num2str(A) num2str(B+1) num2str(CC)];
    Airfoils{1,4} = [num2str(A) num2str(B) num2str(CC+1)];
    
    Airfoils{2,1} = [num2str(A+1) num2str(B) num2str(CC)];
    Airfoils{2,2} = [num2str(A+2) num2str(B) num2str(CC)];
    Airfoils{2,3} = [num2str(A+1) num2str(B+1) num2str(CC)];
    Airfoils{2,4} = [num2str(A+1) num2str(B) num2str(CC+1)];
    
    Airfoils{3,1} = [num2str(A) num2str(B+1) num2str(CC)];
    Airfoils{3,2} = [num2str(A+1) num2str(B+1) num2str(CC)];
    Airfoils{3,3} = [num2str(A) num2str(B+2) num2str(CC)];
    Airfoils{3,4} = [num2str(A) num2str(B+1) num2str(CC+1)];
    
    Airfoils{4,1} = [num2str(A) num2str(B) num2str(CC+1)];
    Airfoils{4,2} = [num2str(A+1) num2str(B) num2str(CC+1)];
    Airfoils{4,3} = [num2str(A) num2str(B+1) num2str(CC+1)];
    Airfoils{4,4} = [num2str(A) num2str(B) num2str(CC+2)];
    
    for r = 1:4
        for c = 1:4
            if length(Airfoils{r,c}) == 3
                Airfoils{r,c} = [Airfoils{r,c}(1:2) '0' Airfoils{r,c}(3)];
            end
        end
    end
end