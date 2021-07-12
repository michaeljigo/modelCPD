function newEcc = convertEcc(expEcc,expDist,newDist)
   % Purpose:  This function will take the eccentricities used in an experiment and the distance
   %           at which the monitor was placed from the subject, and convert those eccentricities
   %           to the expected eccentricities if that same monitor (i.e., same physical size) was 
   %           moved to a different distance from the subject. 
   %
   % Input:    expEcc     --    eccentricities used in the experimental setup
   %           expDist    --    distance from subject to monitor used in experiment
   %           newDist    --    new distance from subject to monitor
   %
   % Output:   newEcc     --    new eccentricities at the new monitor distance 

newEcc = 2*atand((expDist*tand(expEcc./2))./newDist);
