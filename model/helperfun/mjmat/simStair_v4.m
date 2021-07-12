% Fourth iteration of staircase simulations. The goal this time is to make a 
% very general function (or set of functions) such that future simulations will 
% be easier to conduct.

function stair = simStair_v4(subj,stair,simul)

initStair = stair;
for s = 1:simul.n
   stair(s) = initStair;

   % specify field that holds stimulus level
   switch stair(s).type
      case 'downup'
         levelField = 'currentLevel';
      case 'maxlikelihood'
         levelField = 'xCurrent';
   end

   % run staircase
   for t = 1:simul.trials
      thisLevel = stair(s).(levelField);
      % constrain contrast subject sees to be within 0 and 1
      if stair(s).(levelField)<0
         thisLevel = 0;
      elseif stair(s).(levelField)>1
         thisLevel = 1;
      end
      p_correct = evalPF(subj.pf,thisLevel,subj.params);
      correct = p_correct>rand(1);

      % update staircase
      switch stair(s).type
         case 'downup'
            stair(s) = nDown1Up(stair(s),correct);
         case 'maxlikelihood'
            stair(s) = usePalamedesStaircase(stair(s),correct);
      end
   end
end
