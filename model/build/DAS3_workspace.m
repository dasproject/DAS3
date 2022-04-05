function angles = DAS3_workspace(input1,input2,input3)
% Find shoulder workspace for DAS3
% Estimates SCy, SCz, SCx, ACy, ACz, ACx, GHy, GHz and GHyy based on thoracohumeral
% angles, using a regression model of the shoulder rhythm
% 
% 29/03/22: Edit by D Blana to remove elbow angles that are unused

angles=[];

if nargin==1
    das2angles = load(input1);
    allx = das2angles.GHdata;
    for index=1:size(allx,1)
        x = estimate_shoulder_angles(allx(index,1)*pi/180,allx(index,2)*pi/180,allx(index,3)*pi/180);
        angles=[angles x'];  
    end
    angles = angles';

elseif nargin==3
    GHy = input1;
    GHz = input2;
    GHyy = input3;
    for hum_thory=GHy
        for hum_thorz=GHz
            for hum_thoryy=GHyy
                x = estimate_shoulder_angles(hum_thory*pi/180,hum_thorz*pi/180,hum_thoryy*pi/180);
                angles=[angles x'];               
            end
        end
    end
    angles = angles';
end
