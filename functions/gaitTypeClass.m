function gaitType = gaitTypeClass(vel, stRate, gaitClass)
%% Identify gait type using a trained LDA classifier.  This will be more
% robust for shuffle-runners, older adults and speed walkers.  gaitClass
% represents an LDA object which has been trained on 839 test sets of
% walking and running, and validated on ~2000 sets of walking and running.
testSet = [vel stRate];
import classreg.learning.classif.CompactClassificationDiscriminant
%load('gaitClass.mat','gaitClass')


label = predict(gaitClass,testSet);
%label returned as cell
gaitType = label{1};

end