 function [Total_Mass, CMx, CMy, CMz] = Center_Of_Mass_Calculator(inFile)
%     inFile = "C:\TUAS\carbon-copy\MATLAB\CM_Calc\carbon_copy_mass_data.xlsx";
    inSheet = readcell(inFile);
    Total_Mass = 0;
    Total_xTrans = 0;
    Total_yTrans = 0;
    Total_zTrans = 0;
    % inSheet(isnan(inSheet))= 0;

        for i = 2:100
            mass = inSheet{i,5};
            xPosition = inSheet{i,6};
            yPosition = inSheet{i,7};
            zPosition = inSheet{i,8};

                Total_Mass =  Total_Mass + mass;

             xTrans = mass * xPosition;
             Total_xTrans = Total_xTrans + xTrans;
                CMx = Total_xTrans / Total_Mass;

             yTrans = mass * yPosition;
             Total_yTrans = Total_yTrans + yTrans;
                CMy = Total_yTrans / Total_Mass;

             zTrans = mass * zPosition;
             Total_zTrans = Total_zTrans + zTrans;
                CMz = Total_zTrans / Total_Mass;
        end

    % T = table(Total_Mass, CMx, CMy, CMz);
    % writetable(T,inFile,'Sheet', 2,'Range','A1');
    M = [Total_Mass, CMx, CMy, CMz];
    writematrix(M,"C:\TUAS\carbon-copy\MATLAB\CM_Calc\COMC_Output.txt","Delimiter","\t")
 end
 
