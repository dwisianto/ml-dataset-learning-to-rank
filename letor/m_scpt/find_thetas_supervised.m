%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcualte the MEAN of Kendall's Tau rank correlation for each feature 
% for all partial ranking lists based on target labels. 
% It assumes higher ranks' feature values > lower ranks' feature values. 
%        corr = (#agree - #disagree)/ #legal pairs
% Note: 
%    1. smaller numbers means higher rank (i.e., 0 >> 1 >> 2)
%    2. ignore missing values (INF)  
%    3. set negative values to 0
%    4. legal pairs means pairs with different ranks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[correlation_means] = find_thetas_supervised(targets, data)
    nqueries = size(data, 1);
    nlists = size(data{1, 1}, 2);

    rank_correlation = zeros(nqueries, nlists);
    npairs = zeros(nqueries, 1);
    % For each partial list, calculate its Kendall's Tau rank correlation
    for i = 1:size(data, 1)
        [rank_correlation(i, :), npairs(i, 1)] = kendall_tau_rank_correlation(targets{i, 1}, data{i, 1}, nlists);
    end%end-for each partial ranking list
    
    correlation_means = zeros(1, nlists);
    % For each feature, get averaged rank correlation
    for i = 1:nlists
        %compute avergate Kendall's tau ignoring missing values (Inf)
        correlation_means(1, i) = mean(rank_correlation(rank_correlation(:, i) ~= Inf, i), 1);
    end%end-for
    % Set negative mean values to 0
    correlation_means(1, correlation_means(1, :) < 0) = 0;
end%end-function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Kendall's Tau rank correlation for each partial ranking list
% based on target labels. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[rank_correlation, npairs] = kendall_tau_rank_correlation(targets, partial_rankings, nlists)

    nagree = zeros(1, nlists);
    ndisagree = zeros(1, nlists);
    npairs = 0;
    ndocs = size(targets, 1);

    % Enumerate all possible pairs
    for i = 1:ndocs - 1
        for j = i+1:ndocs
            % if rank is the same, skip this round
            if targets(i, 1) == targets(j, 1)
                continue;
            end%end-if
            % get # pairs such that their ranks are different (legal pairs)
            npairs = npairs + 1;

            % if i's rank is behind the j's rank ==> 
            % all i's feature values < all 2nd rank's feature values
            % Note: if either value is 0, skip
            if targets(i, 1) > targets(j, 1)
                for k = 1:nlists
                    % if either one has missing value (0), skip 
                    if partial_rankings(i, k) == 0 || partial_rankings(j, k) == 0
                        continue;
                    end%end-if
                    % 1st values should be smaller than 2nd values 
                    if partial_rankings(i, k) < partial_rankings(j, k)
                        nagree(1, k) = nagree(1, k) + 1;
                    else
                        ndisagree(1, k) = ndisagree(1, k) + 1;
                    end%end-if-else
                end%end-for
            % if 1st rank precceds the 2nd rank ==> 
            % all 1st rank's feature values > all 2nd rank's feature values
            % Note: if either value is 0, skip
            else
                for k = 1:size(partial_rankings, 2)
                    if partial_rankings(i, k) == 0 || partial_rankings(j, k) == 0
                        continue;
                    end%end-if
                    % 1st values should be greater than 2nd values
                    if partial_rankings(i, k) > partial_rankings(j, k)
                        nagree(1, k) = nagree(1, k) + 1;
                    else
                        ndisagree(1, k) = ndisagree(1, k) + 1;
                    end%end-if-else
                end%end-for
            end%end-if-else
        end%end-for
    end%end-for

    % Initialize the rank correlation lists to INF
    rank_correlation = inf(1, nlists);
    % If no legal pairs present, return 
    if npairs == 0
        return;
    end%end-if

    % For each feature, caluculate its rank correlation
    for i = 1:nlists
        % if either agreement or disagreement list is non-empty
        if nagree(1, i) > 0 || ndisagree(1, i) > 0
            % corr = (#agree - #disagree)/ #legal pairs
            rank_correlation(1, i) = (nagree(1, i) - ndisagree(1, i)) ./ npairs;
        end%end-if
    end%end-for
end%end-function
