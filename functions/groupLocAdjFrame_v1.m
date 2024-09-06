% 181029 Tianben Ding

% Group localizations in adjacent frames 

function [posOld,photonOld,precOld,onTimeOld,groupPos,posNonGroupAll,photonNonGroupAll,precNonGroupAll,idNonGroupAll,onTimeNonGroupAll] = ...
    groupLocAdjFrame_v1(posNew,photonNew,precNew,posOld,photonOld,precOld,onTimeOld,groupPos,...
    posNonGroupAll,photonNonGroupAll,precNonGroupAll,idNonGroupAll,onTimeNonGroupAll,groupCoe,h)

% INPUT
% - localization information in the new frame
% POSNEW: localized emitter position in the most recent frame (or working
% frame), [xpos;ypos], 2*n matrix, where n is an arbitrary number of
% localization in the frame
% PHOTONNEW: photon number of each localization, index corresponding to
% posNew, 1*n array
% PRECNEW: localization precision of each localization, index corresponding
% to precNew, 1*n array
% H: current frame index

% - stored localization information up to the immediately preceding frame
% POSOLD: localized emitter position in the immediately preceding frame, 
% [xpos;ypos], 2*m matrix, where m is an arbitrary number of localization 
% in the frame
% PHOTONOLD: photon number of each localization, index corresponding to
% posOld, 1*m array
% PRECOLD: localization precision of each localization, index corresponding
% to precOld, 1*m array
% ONTIMEOLD: emitter on-time of each localization in the immediately 
% preceding frame, unit in frame, index corresponding to precOld, 1*m array
% GROUPPOS: all grouped localization up to the immediately preceding frame,
% only include grouped localization in adjacent frames, not in sequence,
% for checking individual grouping quality, 4*l matrix, where l represents
% all grouped localization pairs so far
% POSNONGROUPALL: stored localization information that are not grouped with
% other localizations up to the immediately preceding frame, 2*k matrix,
% where k represents all localizations that are not grouped up to this
% point
% PHOTONNONGROUPALL, PRECNONGROUPALL, IDNONGROUPALL, ONTIMENONGROUPALL:
% similar as posNonGroupAll, but include photons, localization precision,
% localization id, and emitter on-time instead. Sizes of these horizontal
% arrays should be matched with each other, 1*k arrays
% GROUPCOE: grouping coefficient, unit in localization precision, define
% the radii of grouping circles, try 2~3

% OUTPUT
% POSOLD, PHOTONOLD, PRECOLD, ONTIMEOLD, GROUPPOS, POSNONGROUPALL,
% PHOTONNONGROUPALL, PRECNONGROUPALL, IDNONGROUPALL, ONTIMENONGROUPALL:
% updated input variables after the grouping procedure


% initialize ungrouped index
ungroupedIndOld = ones(1,length(photonOld));
ungroupedIndNew = ones(1,length(photonNew));

posOldGroup = [];
photonOldGroup = [];
precOldGroup = [];
onTimeOldGroup = [];

if ~isempty(photonOld) && ~isempty(photonNew)
    posOldMeshX = repmat(posOld(1,:),length(posNew(1,:)),1);
    posOldMeshY = repmat(posOld(2,:),length(posNew(2,:)),1);
    posNewMeshX = repmat(posNew(1,:).',1,length(posOld(1,:)));
    posNewMeshY = repmat(posNew(2,:).',1,length(posOld(2,:)));
    % calculate Euclidean distance between all SMs in two successive
    % frames (old and new frames)
    allVectorDis = hypot(posOldMeshX-posNewMeshX,posOldMeshY-posNewMeshY);
    
    precMatrix = max(repmat(precOld,length(precNew),1),repmat(precNew.',1,length(precOld)));
    % pick a larger precision for the following grouping, precMatrix
    % contains precisions of all possible localization combinations in the
    % successive frames
    
    candInd = allVectorDis < (groupCoe*precMatrix);
    
    while sum(sum(candInd)) > 0
        candidate = nan(size(allVectorDis));
        candidate(candInd) = allVectorDis(candInd);
        [miniR,minIndR] = nanmin(candidate,[],2);
        [~,minIndC] = nanmin(miniR);
        pairIndSM = [minIndR(minIndC);minIndC];
        groupPos = [groupPos [posOld(:,pairIndSM(1));posNew(:,pairIndSM(2))]];
        posOld(:,pairIndSM(1)) = posNew(:,pairIndSM(2));
        posOldGroup = [posOldGroup posOld(:,pairIndSM(1))];
        photonOld(pairIndSM(1)) = photonOld(pairIndSM(1)) + photonNew(pairIndSM(2));
        photonOldGroup = [photonOldGroup photonOld(pairIndSM(1))];
        precOld(pairIndSM(1)) = precNew(pairIndSM(2));
        precOldGroup = [precOldGroup precOld(pairIndSM(1))];
        onTimeOld(pairIndSM(1)) = onTimeOld(pairIndSM(1)) + 1;
        onTimeOldGroup = [onTimeOldGroup onTimeOld(pairIndSM(1))];
        ungroupedIndOld(pairIndSM(1)) = 0;
        ungroupedIndNew(pairIndSM(2)) = 0;
        candInd(minIndC,:) = 0;
        candInd(:,minIndR(minIndC)) = 0;
    end
end

% store information of localizations that are not grouped with localizations in the new frame
posNonGroupAll = [posNonGroupAll posOld(:,logical(ungroupedIndOld))];
photonNonGroupAll = [photonNonGroupAll photonOld(logical(ungroupedIndOld))];
precNonGroupAll = [precNonGroupAll precOld(logical(ungroupedIndOld))];
idNonGroupAll = [idNonGroupAll (h-1)*ones(1,sum(ungroupedIndOld))];
onTimeNonGroupAll = [onTimeNonGroupAll onTimeOld(logical(ungroupedIndOld))];

% replace "old" frame information with grouped localizations and the new frame information
posOld = [posOldGroup posNew(:,logical(ungroupedIndNew))];
photonOld = [photonOldGroup photonNew(logical(ungroupedIndNew))];
precOld = [precOldGroup precNew(logical(ungroupedIndNew))];
onTimeOld = [onTimeOldGroup ones(1,sum(ungroupedIndNew))];
end