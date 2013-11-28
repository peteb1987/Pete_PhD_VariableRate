function [ ind ] = most_recent_changepoint( list, time )
%MOST_RECENT_CHANGEPOINT Find the index of the most recent changepoint in a
%list.

ind = find(list==max(list(list<=time)),1);

end

